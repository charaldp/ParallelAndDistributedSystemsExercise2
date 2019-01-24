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
    int i,j,k,l,team,processId,noProcesses,*child_Id,*child_num,size,partLength,coordSize,vantagePoint,
    *tempVantagePoints,*localTreeVPs,*commonTreeVPs,count,multiplicity,indexOffset;
    float median,*tempMedians,*distances,*localTreeMedians,*commonTreeMedians;
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
        //The extra coordinate (+ 1) is used to keep the iterations on which this point was chosen as vantage point !!!
        pointsCoords[i] = (floatType *)calloc((coordSize + 1), sizeof(floatType));
        //calloc() is used in order to set the initial iteration bits chosen as vantage point to zero(explained at lines 104 & 238,241)
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
        generatePoints(pointsCoords,partLength,coordSize ,processId);
    }
    
    //Reading data from a file using rowOffset causes asynchrony
    MPI_Barrier(MPI_COMM_WORLD);

    int l_parallel_max = log(noProcesses)/log(2);
    int l_max = log(size)/log(2);
    childComm = (MPI_Comm*)malloc( l_parallel_max * sizeof(MPI_Comm));
    child_num = (int*)malloc( l_parallel_max * sizeof(int));
    child_Id = (int*)malloc( l_parallel_max * sizeof(int));
    commonTreeMedians = (float*)malloc(((1<<l_parallel_max) - 1)*sizeof(float));
    if(processId == 0)
        tempMedians = (float*)malloc(noProcesses * sizeof(float));

    team = processId/noProcesses;

    for(l = 0,indexOffset = 0,childComm[l] = MPI_COMM_WORLD,child_Id[l] = processId,child_num[l] = noProcesses;
        l < l_parallel_max;
        indexOffset += 1<<l,l++){
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
            // Add a bit for each iteration l the chosen point was the vantage point in the extra coordinate
            pointsCoords[vantagePoint][coordSize] += 1<<l;
            #ifdef DEBUG
            printf("Iteration %d Team%dVantagePoint = %d\n",l,team,vantagePoint);
            #endif
            for(i = 0;i<coordSize;i++){
                vantagePointCoords[i] = pointsCoords[vantagePoint][i];
            }
            
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
            median = masterPart(child_num[l],child_Id[l],size/(1<<l),partLength,distances,childComm[l]) + 1/1.0e4;
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
        MPI_Gather(&median, 1, MPI_FLOAT, tempMedians, 1, MPI_FLOAT, 0,MPI_COMM_WORLD);
        
        for(j = 0;j < 1<<l;j++){
            if( processId == 0 ){
                commonTreeMedians[indexOffset + j] = tempMedians[j * (noProcesses/(1<<l))];
            }

            MPI_Bcast(&commonTreeMedians[indexOffset + j],1,MPI_FLOAT,0 ,MPI_COMM_WORLD);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        calculateDistances(distances,pointsCoords,vantagePointCoords,partLength,coordSize);

        // Count points that need exchange
        count = 0;
        
        if(child_Id[l] >= child_num[l]/2){
            for(i = 0;i < partLength;i++)
                if(distances[i] <= median)
                    count++;
        }else{
            for(i = 0;i < partLength;i++)
                if(distances[i] > median)
                    count++;
        }


        #ifdef DEBUG_MAIN
            printf("l = %d,count[%d] = %d\n",l,child_Id[l],count);
        #endif
        MPI_Barrier(MPI_COMM_WORLD);
        transferPoints(distances,median,pointsCoords,partLength,coordSize + 1,child_Id[l],child_num[l],childComm[l],count);
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

    localTreeMedians = (float*)malloc(((1<<(l_max - l_parallel_max)) - 1) * sizeof(float));
    localTreeVPs = (int*)malloc(((1<<(l_max - l_parallel_max)) - 1) * sizeof(int));
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

        for(i = 0;i < multiplicity;i++){
            srand(time(NULL)*(l+1)*(processId+1)*(i+1));
            tempVantagePoints[i] = i * (partLength/multiplicity) + (rand() % (partLength/(multiplicity*2)));
            pointsCoords[tempVantagePoints[i]][coordSize] += 1<<l;
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
            localTreeMedians[indexOffset + i] = tempMedians[i];
            localTreeVPs[indexOffset + i] = tempVantagePoints[i];
        }
        
    }

    gettimeofday(&first, &tzp);
    int index,*allIndices;
    commonTreeVPs = (int*)malloc(((1<<l_parallel_max) - 1)*sizeof(int));
    if(processId == 0)
        allIndices = (int*)malloc(noProcesses * sizeof(int));
    if(processId == 0)
    printf("vpIndex = |");
    for(l = 0,indexOffset = 0;l < l_parallel_max;l++,indexOffset += i){
        k = -1;
        //Create the mask i in order to obtain the vantage points of the current iteration
        i = 1<<l;
        for(j = 0;j < partLength;j++)
            if(((int)(pointsCoords[j][coordSize]) & i) > 0) // The appliance of this mask satisfies this condition only for vantage points of this iteration
                k = processId*partLength + j;

        MPI_Gather(&k, 1, MPI_INT, allIndices, 1, MPI_INT, 0,MPI_COMM_WORLD);
        if(processId == 0)
            index = 0;
        for(j = 0;j < i;j++){
            if(processId == 0){
                while(allIndices[index] == -1){index++;}
                commonTreeVPs[indexOffset + j] = allIndices[index];
                index++;
            }
            MPI_Bcast(&commonTreeVPs[indexOffset + j],1,MPI_INT,0,MPI_COMM_WORLD);
            if(processId == 0)
                printf("%d|",commonTreeVPs[indexOffset + j]);
        }
        if(processId == 0)
            printf("\n");
    }
    gettimeofday(&second, &tzp);
    vpTimeSum += (double)((second.tv_usec - first.tv_usec)/1.0e6
                + second.tv_sec - first.tv_sec);
    MPI_Barrier(MPI_COMM_WORLD);
    printf("Total Time elapsed for creation of vptree on process %d is %f sec\n",processId,vpTimeSum);

    if(processId == 0){
        printf("medians = |");
        for(l = 0,indexOffset = 0;l < l_parallel_max;indexOffset += 1 << l,l++){
            for(j = 0;j < 1<<l;j++)
                printf("%f|",commonTreeMedians[indexOffset + j]);
            printf("\n");
        }
    }

    int co_op,co_opFlag,procCommIt,*co_opFlags,process,neighboursFound,**allPointsKneighbours,interestPoint,partSize,totalNeighboursFound;
    floatType *interestPointCoords,**neighboursCoords;
    neighboursCoords = (floatType**)malloc(sizeof(floatType*));
    interestPointCoords = (floatType*)malloc(coordSize * sizeof(floatType));
    co_opFlags = (int*)malloc(noProcesses * sizeof(int));
    // All kNN-Search using Vantage Point Tree structure
    for(k = 2;k <= 256;k = k << 1){
        allPointsKneighbours = (int**)malloc(size * sizeof((int*)malloc(k * sizeof(int))));
        for(i = 0;i < size;i++)
            allPointsKneighbours[i] = (int*)malloc(k * sizeof(int));

        for(process = 0;process < noProcesses;process++){
            if(process == processId){
                for(interestPoint = 0;interestPoint < partLength;interestPoint++){
                    totalNeighboursFound = 0;
                    interestPointCoords = pointsCoords[interestPoint];
                    for(l = l_max-1,indexOffset = (1 << (l_max - l_parallel_max)) - 1;l >= l_parallel_max;l--){
                        multiplicity = 1<<(l - l_parallel_max);
                        partSize = partLength/multiplicity;
                        for(i = (interestPoint/partSize) * partSize;
                            i < ((interestPoint/partSize) + 1) * partSize;
                            i++){
                            if(i != interestPoint){
                                calculateDistances(distances,pointsCoords,interestPointCoords,partSize,coordSize);
                                if(distance(pointsCoords[i],interestPointCoords,coordSize) <=   //If the interest point's distance from the current neighbour is less than the distance from vptree's current hypersphere's bound
                                 abs(localTreeMedians[indexOffset + i/partSize] - distance(pointsCoords[i],pointsCoords[localTreeVPs[indexOffset + i/partSize]%partLength],coordSize)) ){
                                    allPointsKneighbours[totalNeighboursFound] = partLength * processId + i;
                                    totalNeighboursFound++; // A certain neighbour has been found
                                }
                            }
                        }
                        for(procCommIt = 0;procCommIt < noProcesses;procCommIt++)
                            if(procCommIt != processId)
                                MPI_Send(&co_opFlags[procCommIt],1,MPI_INT,procCommIt,1,MPI_COMM_WORLD);

                        MPI_Send(interestPointCoords,coordSize,MPI_floatType,co_op,1,MPI_COMM_WORLD);
                        MPI_Recv(&neighboursFound,1,MPI_INT,co_op,1,MPI_COMM_WORLD,&Stat);
                        if(neighboursFound > 0){
                            neighboursCoords = (floatType**)realloc(neighboursCoords,neighboursFound * coordSize * sizeof(floatType*));

                            MPI_Recv(neighboursCoords,neighboursFound,MPI_floatType,process,1,MPI_COMM_WORLD,&Stat);
                        }

                    }

                }
                
            }else{
                MPI_Recv(&co_opFlag,1,MPI_INT,process,1,MPI_COMM_WORLD,&Stat);
                if(co_opFlag){
                    MPI_Recv(interestPointCoords,coordSize,MPI_floatType,process,1,MPI_COMM_WORLD,&Stat);
                    //TODO Search on my subtree
                    neighboursCoords = (floatType**)realloc(neighboursCoords,neighboursFound * coordSize * sizeof(floatType*));
                    MPI_Send(&neighboursFound,1,MPI_INT,process,1,MPI_COMM_WORLD);
                    if(neighboursFound > 0)
                        MPI_Send(neighboursCoords,neighboursFound,MPI_floatType,process,1,MPI_COMM_WORLD);
                }
            }
        }
    }

    MPI_Finalize();
    exit(0);
}