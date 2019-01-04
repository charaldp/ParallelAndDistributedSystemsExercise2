/**********************************************************************
 *
 * vptree.c
 * Vantage point tree implementation
 * Charalampos Papadiakos <charaldp@auth.gr>
 * Time-stamp: <2019-01-4>
 *
 **********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <mpi.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "mpiFindMedian.h"

int main(int argc, char **argv)
{
    int team,processId,noProcesses,*child_Id,*child_num,size,partLength,i,l,coordSize,vantagePoint,*counts;
    float median,*distances,*vpMedianDistances;
    floatType **pointsCoords,*vantagePointCoords;
    MPI_Comm *childComm;
    char *dataset;
    MPI_Init (&argc, &argv);                        // starts MPI //
    MPI_Comm_rank (MPI_COMM_WORLD, &processId);     // get current process id //
    MPI_Comm_size (MPI_COMM_WORLD, &noProcesses);   // get number of processes //
    size = 1 << atoi(argv[1]);
    coordSize = atoi(argv[2]);
    
    if (coordSize < 1 || argc < 3 || argc > 4){
        if(processId==0)
            printf("Usage:mpiexec/mpirun -np/-n Num_Threads(Power of 2) %s <log_2(Size)> <Dimensions> <Dataset_Filename>(Optional)\n",argv[0]);
        MPI_Finalize();
        return 0;
    }else if(argc == 4){
        dataset = argv[3];
    }
    counts = (int *)malloc(3 * sizeof(int));
    partLength = size/noProcesses;
    assert(size%noProcesses == 0);
    //Allocating memory for data of each process
    pointsCoords = (floatType **)malloc(partLength * sizeof(floatType*));
    for(i = 0; i < partLength; i++){
        pointsCoords[i] = (floatType *)malloc((coordSize-1) * sizeof(floatType));
    }
    distances = (float *)malloc(partLength * sizeof(float));
    vantagePointCoords = (floatType*)malloc(coordSize * sizeof(floatType));

    // HIGGS has 11.000.000x29 cols rows and CUSY.csv has 500.000x30 cols
    // cities.txt requires an column offset of 4 and has 2 useful cols
    if(argc == 4){
        printf("Dataset : %s\n",dataset);
        if(strcmp(dataset,"cities.txt")==0)
            read_csv(partLength, coordSize, dataset, pointsCoords,processId * partLength,4);
        else
            read_csv(partLength, coordSize, dataset, pointsCoords,processId * partLength,0);
    }else{
        generatePoints(pointsCoords,partLength,coordSize,processId);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    /*dataPartition*/
    childComm = (MPI_Comm*)malloc( ( log(noProcesses)/log(2) + 1) * sizeof(MPI_Comm));
    child_num = (int*)malloc( ( log(noProcesses)/log(2) + 1) * sizeof(int));
    child_Id = (int*)malloc( ( log(noProcesses)/log(2) + 1) * sizeof(int));
    MPI_Comm Current_Comm = MPI_COMM_WORLD;
    for(l = 0;l < log(noProcesses)/log(2); l++){

        team = processId/(noProcesses/(1<<l));
        
        MPI_Comm_split(Current_Comm, team, processId, &childComm[l]);
        MPI_Comm_rank(childComm[l], &child_Id[l]);     // get current process id //
        MPI_Comm_size(childComm[l], &child_num[l]); 
        printf("l = %d, I %d from Team - %d, was split to %d of %d children\n",l,processId,team,child_Id[l],child_num[l]);

        if(child_Id[l]==0)
        {
            srand(time(NULL)*(team+1)*(l+1)*(processId+1));
            vantagePoint = rand() % partLength;
            printf("Iteration %d Team%dVantagePoint = %d\n",l,team,vantagePoint);
            for (i = 0; i < coordSize; i++) {
                vantagePointCoords = pointsCoords[vantagePoint];
            }

                //MPI_Sendnudes(&partLength,1,MPI_INT,processId,MPI_COMM_WORLD); /*Test ( ͡° ͜ʖ ͡°)*/
                
                MPI_Bcast(vantagePointCoords,coordSize,MPI_floatType,child_Id[l],childComm[l]);
                calculateDistances(distances,pointsCoords,vantagePointCoords,partLength,coordSize);
        }
        else
        {

            MPI_Bcast(vantagePointCoords,coordSize,MPI_floatType,0,childComm[l]);
            
            calculateDistances(distances,pointsCoords,vantagePointCoords,partLength,coordSize);
        }
        #ifdef DEBUG
        for(i=0; i<partLength; i++)
            printf("distBefore %d = |%f|\n",i+child_Id[l]*partLength,distances[i]);
        #endif
        if(child_Id[l]==0)
        {
            median = masterPart(child_num[l],child_Id[l],size/(1<<l),partLength,distances,childComm[l],counts);
            printf("l = %d,Median from team %d is: %f\n",l,team,median);
        }
        else{
            slavePart(child_Id[l],partLength,distances,size/(1<<l),child_num[l],childComm[l],counts);
        }
        MPI_Bcast(&median,1,MPI_FLOAT,0,childComm[l]);
        MPI_Barrier(MPI_COMM_WORLD);
        //transferPointsST(distances,median,pointsCoords,partLength,coordSize);
        calculateDistances(distances,pointsCoords,vantagePointCoords,partLength,coordSize);
 
        transferPoints(distances,median,pointsCoords,partLength,coordSize,child_Id[l],child_num[l],childComm[l],counts);
        MPI_Barrier(MPI_COMM_WORLD);
        calculateDistances(distances,pointsCoords,vantagePointCoords,partLength,coordSize);
        #ifdef DEBUG
        for(i=0; i<partLength; i++)
            printf("distAfter %d = |%f|\n",i+processId*partLength,distances[i]);
        #endif
        validationST(median,partLength,distances,counts,processId);
        assert((child_Id[l] < child_num[l]/2 && (counts[1] + counts[2]) == partLength) || (child_Id[l] >= child_num[l]/2 && counts[0] == partLength));
        MPI_Barrier(MPI_COMM_WORLD);
        Current_Comm = childComm[l];
        printf("l = %d,team = %d, Pid %d CounterMax = %d\n",l,team,processId,counts[0]);
    }
    // Go serial
    for(;l<log(size)/log(2);l++){

        calculateDistances(distances,pointsCoords,vantagePointCoords,partLength,coordSize);
        struct timeval first, second, lapsed;
        struct timezone tzp;
        gettimeofday(&first, &tzp);
        median = selection(distances,partLength);
        gettimeofday(&second, &tzp);
        if(first.tv_usec>second.tv_usec)
        {
            second.tv_usec += 1000000;
            second.tv_sec--;
        }
        lapsed.tv_usec = second.tv_usec - first.tv_usec;
        lapsed.tv_sec = second.tv_sec - first.tv_sec;
        validationST(median,partLength,distances,&counts[0],processId);
        printf("Time elapsed: %lu, %lu s\n", lapsed.tv_sec, lapsed.tv_usec);
        printf("l = %d, Median at Serial from Team %d: %f\n",l,team,median);
        printf("Pid %d CounterMax = %d\n",processId,counts[0]);
        //MPI_Finalize();
        //exit(0);
    }
    // knn Search


    MPI_Finalize();
    //free(distances);
    exit(0);

}