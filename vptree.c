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



double vpTimeSum = 0, knnTimeSum = 0;

int main(int argc, char **argv)
{
    int i,k,l,team,processId,noProcesses,*child_Id,*child_num,size,partLength,coordSize,vantagePoint,*tempVantagePoints,*allVantagePoints,*commonTreeVPs,count;
    float median,*tempMedians,*allMedians,*distances,*commonTreeMedians;
    floatType **pointsCoords,*vantagePointCoords;
    MPI_Comm *childComm;
    char *dataset;
    struct timeval first, second;
    struct timezone tzp;
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
        if(processId == 0)
            printf("Reading from Dataset : %s\n",dataset);
        read_csv(partLength, coordSize, dataset, pointsCoords,processId * partLength,
                 strcmp(dataset,"cities.txt") ? 0 : 4);      //Am I reading from cities.txt?
    }else{
        generatePoints(pointsCoords,partLength,coordSize,processId);
    }

    //Reading data from a file using rowOffset causes asynchrony
    MPI_Barrier(MPI_COMM_WORLD);

    int l_parallel_max = log(noProcesses)/log(2);
    childComm = (MPI_Comm*)malloc( l_parallel_max * sizeof(MPI_Comm));
    child_num = (int*)malloc( l_parallel_max * sizeof(int));
    child_Id = (int*)malloc( l_parallel_max * sizeof(int));
    
    team = processId/noProcesses;

    for(l = 0,childComm[l] = MPI_COMM_WORLD,child_Id[l] = processId,child_num[l] = noProcesses;
        l < l_parallel_max;
        l++){
        gettimeofday(&first, &tzp);
        if(l > 0){
            team = processId/(noProcesses/(1<<l));
            MPI_Comm_split(childComm[l-1], team, processId, &childComm[l]);
            MPI_Comm_rank(childComm[l], &child_Id[l]);     // get current process id //
            MPI_Comm_size(childComm[l], &child_num[l]);
        }

        #ifdef DEBUG_MAIN
            printf("l = %d, I %d from Team - %d, was split to %d of %d children\n",l,processId,team,child_Id[l],child_num[l]);
        #endif

        if(child_Id[l]==0)
        {
            srand(time(NULL)*(team+1)*(l+1)*(processId+1));
            vantagePoint = rand() % partLength;
            #ifdef DEBUG
            printf("Iteration %d Team%dVantagePoint = %d\n",l,team,vantagePoint);
            #endif
            vantagePointCoords = pointsCoords[vantagePoint];

            //MPI_Sendnudes(&partLength,1,MPI_INT,processId,MPI_COMM_WORLD);
            MPI_Bcast(vantagePointCoords,coordSize,MPI_floatType,child_Id[l],childComm[l]);
            calculateDistances(distances,pointsCoords,vantagePointCoords,partLength,coordSize);
        }
        else
        {

            MPI_Bcast(vantagePointCoords,coordSize,MPI_floatType,0,childComm[l]);
            calculateDistances(distances,pointsCoords,vantagePointCoords,partLength,coordSize);
        }
        gettimeofday(&second, &tzp);
        vpTimeSum += (double)((second.tv_usec - first.tv_usec)/1.0e6
                    + second.tv_sec - first.tv_sec);
        #ifdef DEBUG_DIST
        for(i = 0; i<partLength; i++)
            printf("distBefore %d = |%f|\n",i+child_Id[l]*partLength,distances[i]);
        #endif
        if(child_Id[l]==0)
        {
            median = masterPart(child_num[l],child_Id[l],size/(1<<l),partLength,distances,childComm[l]);
            #ifdef DEBUG
            printf("l = %d,Median from team %d is: %f\n",l,team,median);
            #endif
        }
        else{
            slavePart(child_Id[l],partLength,distances,size/(1<<l),child_num[l],childComm[l]);
        }
        gettimeofday(&first, &tzp);
        MPI_Bcast(&median,1,MPI_FLOAT,0,childComm[l]);
        MPI_Barrier(MPI_COMM_WORLD);
        //transferPointsST(distances,median,pointsCoords,partLength,coordSize);
        calculateDistances(distances,pointsCoords,vantagePointCoords,partLength,coordSize);

        // Count points that need exchange
        count = 0;
        for(i = 0;i < partLength;i++)
            if((child_Id[l] >= child_num[l]/2 && distances[i] <= median) || (child_Id[l] < child_num[l]/2 && distances[i] > median))
                count++;

        transferPoints(distances,median,pointsCoords,partLength,coordSize,child_Id[l],child_num[l],childComm[l],count);
        gettimeofday(&second, &tzp);
        vpTimeSum += (double)((second.tv_usec - first.tv_usec)/1.0e6
                    + second.tv_sec - first.tv_sec);
        calculateDistances(distances,pointsCoords,vantagePointCoords,partLength,coordSize);
        #ifdef DEBUG_DIST
        for(i=0; i<partLength; i++)
            printf("distAfter %d = |%f|\n",i+processId*partLength,distances[i]);
        MPI_Barrier(MPI_COMM_WORLD);
        #endif
        validationPartition(median,partLength,distances,child_Id[l],child_num[l]);
        MPI_Barrier(MPI_COMM_WORLD);
        #ifdef DEBUG
        printf("l = %d,team = %d, Pid %d count = %d\n",l,team,processId,count);
        #endif
    }
    // Single Thread Operations
    tempMedians = (float*)malloc(sizeof(float));
    tempVantagePoints = (int*)malloc(sizeof(int));
    allMedians = (float*)malloc(sizeof(float));
    allVantagePoints = (int*)malloc(sizeof(int));
    int multiplicity,indexOffset;
    for(l = l_parallel_max, indexOffset = 0;
        l<log(size)/log(2);
        l++, indexOffset += multiplicity){
        gettimeofday(&first, &tzp);
        multiplicity = 1<<(l - l_parallel_max);
        #ifdef DEBUG_MAIN
        printf("multiplicity = %d\n",multiplicity);
        #endif
        tempMedians = (float*)realloc(tempMedians,multiplicity * sizeof(float));
        tempVantagePoints = (int*)realloc(tempVantagePoints,multiplicity * sizeof(int));
        allMedians = (float*)realloc(allMedians,(indexOffset + multiplicity) * sizeof(float));
        allVantagePoints = (int*)realloc(allVantagePoints,(indexOffset + multiplicity) * sizeof(int));

        for(i = 0;i < multiplicity;i++){
            srand(time(NULL)*(l+1)*(processId+1)*(i+1));
            tempVantagePoints[i] = i * (partLength/multiplicity) + (rand() % (partLength/(multiplicity*2)));
        }

        calculateDistancesST(distances,pointsCoords,tempVantagePoints,partLength,multiplicity,coordSize);
        tempMedians = multiSelection(distances, partLength, multiplicity);
        calculateDistancesST(distances,pointsCoords,tempVantagePoints,partLength,multiplicity,coordSize);
        transferPointsST(distances,tempMedians,pointsCoords,partLength,multiplicity);
        gettimeofday(&second, &tzp);
        vpTimeSum += (double)((second.tv_usec - first.tv_usec)/1.0e6
                    + second.tv_sec - first.tv_sec);
        #ifdef DEBUG_VPST
        for(i = 0;i < multiplicity;i++)
            printf("l = %d,Pid = %d, tempMedians[%d] = %f\n",l,processId,i,tempMedians[i]);
        #endif
        calculateDistancesST(distances,pointsCoords,tempVantagePoints,partLength,multiplicity,coordSize);
        validationST(tempMedians,partLength,distances,processId,multiplicity);
        calculateDistancesST(distances,pointsCoords,tempVantagePoints,partLength,multiplicity,coordSize);
        validationPartitionST(tempMedians,partLength,distances,multiplicity);
        for(i = 0;i < multiplicity;i++){
            allMedians[indexOffset + i] = tempMedians[i];
            allVantagePoints[indexOffset + i] = tempVantagePoints[i];
        }
        //MPI_Finalize();
        //exit(0);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    printf("Total Time elapsed for creation of vptree on process %d is %f sec\n",processId,vpTimeSum);
    if(processId == 0)
        printf("Serial VP tree array size: %d\n",indexOffset);
    // All kNN-Search using Vantage Point Tree structure
    for(i = 1;i <= 8;i++){
        k = 1 << i;

    }

    MPI_Finalize();
    //free(distances);
    exit(0);
}