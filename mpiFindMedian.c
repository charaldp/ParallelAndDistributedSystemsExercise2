/*
The MIT License (MIT)

Copyright (c) 2014

Athanassios Kintsakis
Contact
athanassios.kintsakis@gmail.com
akintsakis@issel.ee.auth.gr


Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/
/**********************************************************************
 *
 * mpiFindMedian.c
 * Find median and vantage point tree functions 
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

/***Kills processes that have no values left in their arrays****/
void removeElement(int *array, int *size, int element)
{
    int i;
    int flag=0;
    for(i=0;i<*size;i++)
    {
        if(flag==1)
            array[i]=array[i+1];
        if(array[i]==element&& flag==0)
        {
            array[i]=array[i+1];
            flag=1;
        }
    }
    *size=*size-1;
}

/****Swaps two values in an array****/
void swap_values(float *array,int x,int y){
    float temp;
    temp=array[x];
    array[x]=array[y];
    array[y]=temp;
}

/*** Generate a random point set for each process ***/
void generatePoints(floatType **pointsCoords,int pointsLength,int coordSize,int cal){
    srand((cal+1)*time(NULL));     //Generate number to fill the array
    int i,j;
    for(i = 0;i < pointsLength;i++){
        //printf("point %d = [ ",i);
        for(j = 0;j < coordSize;j++){
            pointsCoords[i][j] = (floatType)(rand()) / (floatType)(RAND_MAX) * (floatType)(5000.0)-2500.0;//(float)(300+i)*(float)(5+i)*(float)(i-3)/(float)(partLength)/7.0;
            //printf("%f,",pointCoords[i][j]);
        }
    //printf("\b]\n");
    }
}

/****Calculate distance****/
void calculateDistances(float *distances,floatType **pointCoords,floatType *vantagePointCoords,int pointsLength,int cordSize){
    int i,j;
    for(i=0;i<pointsLength;i++){
        distances[i] = 0;
        for(j=0;j<cordSize;j++){
            distances[i] = distances[i] + (float)pow((double)pointCoords[i][j] - vantagePointCoords[j],2.0);
        }
        distances[i] = sqrt((double)distances[i]);
    }

}

/***Validates the stability of the operation****/
void validation(float median,int partLength,int size,float *numberPart,int processId,MPI_Comm Current_Comm,int *counts)
{   
    MPI_Bcast(&median,1,MPI_FLOAT,0,Current_Comm);
	int countMin=0;
    int countMax=0;
    int countEq=0;
    int sumMax,sumMin,sumEq,i;
    for(i=0;i<partLength;i++)
    {
        if(numberPart[i]>median)
            countMax++;
        else if(numberPart[i]<median)
            countMin++;
        else{
            countEq++;
        }

    }
    counts[0] = countMax;
    counts[1] = countEq;
    counts[2] = countMin;
    printf("I am %d, countMax = %d\n",processId,countMax);
    MPI_Reduce(&countMax,&sumMax,1,MPI_INT,MPI_SUM,0,Current_Comm);
    MPI_Reduce(&countMin,&sumMin,1,MPI_INT,MPI_SUM,0,Current_Comm);
    MPI_Reduce(&countEq,&sumEq,1,MPI_INT,MPI_SUM,0,Current_Comm);
    if(processId==0)
    {
        if((sumMax<=size/2)&&(sumMin<=size/2))  //Checks if both the lower and higher values occupy less than 50% of the total array.
            printf("VALIDATION PASSED!\n");
        else
            printf("VALIDATION FAILED!\n");


	printf("Values greater than median: %d\n",sumMax);
    printf("Values equal to median: %d\n",sumEq);
   	printf("Values lower than median: %d\n",sumMin);
    }
}

/***Validates the stability of the operation (Single Threaded)****/
void validationST(float median,int size,float *numberPart,int *counts,int processId)
{
	int countMin=0;
    int countMax=0;
    int countEq=0;
    int i;
    for(i=0;i<size;i++)
    {
        if(numberPart[i]>median)
            countMax++;
        else if(numberPart[i]<median)
            countMin++;
        else
            countEq++;
    }
    if((countMax<=size/2)&&(countMin<=size/2))  //Checks if both the lower and higher values occupy less than 50% of the total array.
        printf("VALIDATION PASSED!\n");
    else
        printf("VALIDATION FAILED!\n");


	printf("Values greater than median %d: %d\n",processId,countMax);
    printf("Values equal to median %d: %d\n",processId,countEq);
    printf("Values lower than median %d: %d\n",processId,countMin);
    counts[0] = countMax;
    counts[1] = countEq;
    counts[2] = countMin;
}

void validationPartition(float median,int size,float *numberPart,int *counts,int processId)
{
    int countMin=0;
    int countMax=0;
    int countEq=0;
    int i;
    for(i=0;i<size;i++)
    {
        if(numberPart[i]>median)
            countMax++;
        else if(numberPart[i]<median)
            countMin++;
        else
            countEq++;
    }
    if((countMax<=size/2)&&(countMin<=size/2))  //Checks if both the lower and higher values occupy less than 50% of the total array.
        printf("VALIDATION PASSED!\n");
    else
        printf("VALIDATION FAILED!\n");

    printf("Values greater than median %d: %d\n",processId,countMax);
    printf("Values equal to median %d: %d\n",processId,countEq);
    printf("Values lower than median %d: %d\n",processId,countMin);
    counts[0] = countMax;
    counts[1] = countEq;
    counts[2] = countMin;
}

/****Part executed only by the Master Node****/
float masterPart(int noProcesses,int processId,int size,int partLength,float *numberPart,MPI_Comm Current_Comm,int *countOuter) //MASTER NODE CODE
{
    int elements,i,keepBigSet,sumSets,finalize,randomNode,k,*activeNodes;
    int endSmall=0;
    int dropoutFlag=0;
    int endBig=0;
    float *arraySmall,*arrayBig,*arrayToUse,median,pivot,tempPivot;
    int activeSize=noProcesses;
    int stillActive=1;
    int oldSumSets=-1;
    int checkIdentical=0;
    int useNewPivot=0;
    int *pivotArray;
    k=(int)size/2+1; //It is done so in order to find the right median in an even numbered array.
    elements=partLength;
    activeNodes=(int *)malloc(noProcesses*sizeof(int));  //we create the array that contains the active nodes.
    arrayToUse=numberPart;
    pivotArray=(int *)malloc(noProcesses*sizeof(int));  //Used for special occasions to gather values different than the pivot.
    for(i=0;i<activeSize;i++)
    {
        activeNodes[i]=i;
    }
    int randomCounter=0;
    int randomCounter2=0;
    struct timeval first, second, lapsed;
    struct timezone tzp;
    gettimeofday(&first, &tzp);
    for(;;)   //Begin the infinite loop until the median is found.
    {   
        int counter=0;
        useNewPivot=0;
        if(stillActive==1&&checkIdentical!=0)  //If i still have values in my array and the Sumed Big Set is identical to the previous one, check for identical values.
        {
            for(i=0;i<elements;i++)
            {
                if(pivot==arrayToUse[i])
                    counter++;
                else
                {
                    useNewPivot=1;
                    tempPivot=arrayToUse[i];
                    break;
                }
            }
        }
        if(checkIdentical!=0)
        {
            int useNewPivotMax=0;
	        MPI_Reduce(&useNewPivot,&useNewPivotMax,1,MPI_INT,MPI_MAX,0,Current_Comm); //FIRST(OPTIONAL) REDUCE : MAX useNewPivot
            if(useNewPivotMax!=1)    //That means that the only values left are equal to the pivot!
            {
                median=pivot;
                finalize=1;
                MPI_Bcast(&finalize,1,MPI_INT,processId,Current_Comm); //FIRST(OPTIONAL) BROADCAST : WAIT FOR FINALIZE COMMAND OR NOT
                gettimeofday(&second, &tzp);
                if(first.tv_usec>second.tv_usec)
                {
                    second.tv_usec += 1000000;
                    second.tv_sec--;
                }
                lapsed.tv_usec = second.tv_usec - first.tv_usec;
                lapsed.tv_sec = second.tv_sec - first.tv_sec;
                printf("Time elapsed: %lu, %lu s\n", lapsed.tv_sec, lapsed.tv_usec);
                validation(median,partLength,size,numberPart,0,Current_Comm,countOuter);
                MPI_Barrier(Current_Comm);
                free(pivotArray);
                return median;
            }
            else
            {
                finalize=0;
                int useit=0;
                randomCounter2++;
                MPI_Bcast(&finalize,1,MPI_INT,0,Current_Comm);
                MPI_Gather(&useNewPivot, 1, MPI_INT, pivotArray, 1, MPI_INT, 0, Current_Comm); //Gather every value and chose a node to change the pivot.
                for(i=0;i<activeSize;i++)
                {
                    if(pivotArray[i]==1)
                    {
                        if((randomCounter2>1)&&(randomNode!=activeNodes[i]))  //Check if the same node has already been used in a similar operation.
                        {
                            useit=1;
                            randomNode=activeNodes[i];
                            randomCounter2=0;
                            break;
                        }
                        else if(randomCounter2<2)
                        {
                            useit=1;
                            randomNode=activeNodes[i];
                            break;
                        }
                    }
                }
                if(useit!=0)
                    useNewPivot=1;
                else
                    useNewPivot=0;
            }
        }
        if(useNewPivot!=0)
            MPI_Bcast(&randomNode,1,MPI_INT,0,Current_Comm);  //THIRD(OPTIONAL) BROADCAST : BROADCAST THE SPECIAL NODE
        if(useNewPivot==0)  //if we didnt choose a special Node, choose the node that will pick the pivot in a clockwise manner. Only selects one of the active nodes.
        {
            if(randomCounter>=activeSize)
                randomCounter=0; //Fail safe
            randomNode=activeNodes[randomCounter];
            randomCounter++;			//Increase the counter
            MPI_Bcast(&randomNode,1,MPI_INT,0,Current_Comm);   //FIRST BROADCAST : SENDING randomnode, who will chose
        }
        if(randomNode==processId)  //If i am to choose the pivot.....
	    {
            if(useNewPivot==0)
            {
                srand(time(NULL));
                pivot=arrayToUse[rand() % elements];
                MPI_Bcast(&pivot,1,MPI_FLOAT,0,Current_Comm); //SECOND BROADCAST : SENDING PIVOT   k ton stelnw sto lao
	        }
            else
            {
                MPI_Bcast(&tempPivot,1,MPI_FLOAT,0,Current_Comm); //SECOND BROADCAST : SENDING PIVOT   k ton stelnw sto lao
                pivot=tempPivot;
            }
        }
        else //If not.. wait for the pivot to be received.
            MPI_Bcast(&pivot,1,MPI_FLOAT,randomNode,Current_Comm);  // SECOND BROADCAST : RECEIVING PIVOT
        if(stillActive==1)  //If i still have values in my array.. proceed
        {
            partition(arrayToUse,elements,pivot,&arraySmall,&arrayBig,&endSmall,&endBig);  //I partition my array  // endsmall=number of elements in small array, it may be 0
            // endbig=number of elements in big array, it may be 0
            //arraysmall = Points to the position of the small array.NULL if the array is empty
            //Same for arraybig
        }
        else  //If i'm not active endBig/endSmall has zero value.
        {
            endBig=0;
            endSmall=0;
        }
        sumSets=0;
	    //We add the bigSet Values to decide if we keep the small or the big array
	    MPI_Reduce(&endBig,&sumSets,1,MPI_INT,MPI_SUM,0,Current_Comm);  //FIRST REDUCE : SUM OF BIG
        MPI_Bcast(&sumSets,1,MPI_INT,0,Current_Comm);
        if(oldSumSets==sumSets)
            checkIdentical=1;
        else
        {
            oldSumSets=sumSets;
            checkIdentical=0;
        }
	    //hmetabliti keepBigSet 0 h 1 einai boolean k me autin enimerwnw ton lao ti na kratisei to bigset h to smallset
	    if(sumSets>k)   //an to sumofbigsets > k tote krataw to big SET
	    {
            keepBigSet=1; //to dilwnw auto gt meta tha to steilw se olous
            if(endBig==0)
                dropoutFlag=1; //wraia.. edw an dw oti to bigset mou einai 0.. alla prepei na kratisw to bigset sikwnw auti ti simaia pou simainei tha ginw inactive ligo pio katw tha to deis
            else
            {
                arrayToUse=arrayBig; //thetw ton neo pinaka na einai o big
                elements=endBig; //thetw arithmo stoixeiwn iso me tou big
            }
	    }
	    else if(sumSets<k) //antistoixa an to sumofbigsets < k tote krataw to small set
	    {
		    keepBigSet=0;
		    k=k-sumSets;
		    if(endSmall==0)
                dropoutFlag=1; //antistoixa koitaw an tha ginw inactive..
		    else
		    {
		    	arrayToUse=arraySmall; //dinw times..
		    	elements=endSmall;
		    }
	    }
	    else  //edw simainei k=sumofbigsetes ara briskw pivot k telos
	    {
		    median=pivot;
		    finalize=1; //dilwnw finalaize =1
		    MPI_Bcast(&finalize,1,MPI_INT,0,Current_Comm); //to stelnw se olous, oi opoioi an laboun finalize =1 tote kaloun MPI finalize k telos
		    gettimeofday(&second, &tzp);
            if(first.tv_usec>second.tv_usec)
            {
                second.tv_usec += 1000000;
                second.tv_sec--;
            }
            lapsed.tv_usec = second.tv_usec - first.tv_usec;
            lapsed.tv_sec = second.tv_sec - first.tv_sec;
            printf("Time elapsed: %lu, %lu s\n", lapsed.tv_sec, lapsed.tv_usec);
		    validation(median,partLength,size,numberPart,processId,Current_Comm,countOuter);
            MPI_Barrier(Current_Comm);
            free(pivotArray);
            return median;
        }
        finalize=0; //an den exw mpei sta if den exw steilei timi gia finalize.. oi alloi omws perimenoun na laboun kati, stelnw loipon to 0 pou simainei sunexizoume
        MPI_Bcast(&finalize,1,MPI_INT,0,Current_Comm);	//SECOND BROADCAST : WAIT FOR FINALIZE COMMAND OR NOT
        //edw tous stelnw to keepbigset gia na doun ti tha dialeksoun
	    MPI_Bcast(&keepBigSet,1,MPI_INT,0,Current_Comm);    //THIRD BROADCAST: SEND keepBigset boolean
        if(dropoutFlag==1 && stillActive==1) //edw sumfwna me to dropoutflag pou orisame prin an einai 1 kalw tin sinartisi pou me petaei apo ton pinaka. episis koitaw na eimai active gt an me exei idi petaksei se proigoumeni epanalispi tote den xreiazetai na me ksanapetaksei
        {
            stillActive=0;
            removeElement(activeNodes, &activeSize, processId);
        }
        int flag;
        //edw perimenw na akousw apo ton kathena an sunexizei active h oxi.. an oxi ton petaw.. an einai idi inactive apo prin stelnei kati allo (oxi 1)k den ton ksanapetaw
        for(i=0;i<activeSize;i++)
        {
            if(activeNodes[i]!=processId)
            {
                MPI_Recv(&flag,1,MPI_INT,activeNodes[i],1,Current_Comm,&Stat);  //FIRST RECEIVE : RECEIVE active or not
                if(flag==1)
                    removeElement(activeNodes, &activeSize, activeNodes[i]);
            }
        }
    }
}

/***Executed only by Slave nodes!!*****/
void slavePart(int processId,int partLength,float *numberPart,int size,int noProcesses,MPI_Comm Current_Comm,int *countOuter)  //code here is for the cheap slaves :P
{
	int dropoutflag,elements,i,sumSets,finalize,keepBigSet,randomNode;
    int endSmall=0;
    int endBig=0;
    float *arraySmall,*arrayBig,*arrayToUse,tempPivot,pivot;
	arrayToUse=numberPart;
	elements=partLength;
	int stillActive=1;
	int *pivotArray;
    int oldSumSets=-1;
    int checkIdentical=0;
    int useNewPivot;
    for(;;)
	{
        finalize=0;
        int counter=0;
        useNewPivot=0;
        if(stillActive==1&&checkIdentical!=0)  //If i still have values in my array..   If the Sumed Big Set is identical to the previous one, check for identical values.
        {
            for(i=0;i<elements;i++)
            {
                if(pivot==arrayToUse[i])
                    counter++;
                else
                {
                    useNewPivot=1;
                    tempPivot=arrayToUse[i];
                    break;
                }
            }
        }
        if(checkIdentical!=0)
        {
            int useNewPivotMax=0;
            MPI_Reduce(&useNewPivot,&useNewPivotMax,1,MPI_INT,MPI_MAX,0,Current_Comm);
            MPI_Bcast(&finalize,1,MPI_INT,0,Current_Comm);//an o master apo to keepbigset k apo to count apofasisei oti teleiwsame mou stelnei 1, alliws 0 sunexizoume
            if(finalize==1)
            {
                float median=0;
                validation(median,partLength,size,numberPart,processId,Current_Comm,countOuter);
                MPI_Barrier(Current_Comm);
                return ;
            }
            else
            {
                MPI_Gather(&useNewPivot, 1, MPI_INT, pivotArray, 1, MPI_INT, 0, Current_Comm);
            }
        }
        MPI_Bcast(&randomNode,1,MPI_INT,0,Current_Comm); //FIRST BROAD CAST : RECEIVING RANDOM NODE, perimenw na dw poios einaito done
        if(randomNode!=processId){ //means I am not the one to chose pivot.. so I wait to receive the pivot
            MPI_Bcast(&pivot,1,MPI_FLOAT,randomNode,Current_Comm);	//SECOND BROADCAST : RECEIVING PIVOT
        } 
        else if(randomNode==processId) //I am choosing suckers
        {
            if(useNewPivot==0)
            {
                srand(time(NULL));
                pivot=arrayToUse[rand() % elements];
                
                MPI_Bcast(&pivot,1,MPI_FLOAT,processId,Current_Comm); //SECOND BROADCAST : SENDING PIVOT   k ton stelnw sto lao
            }
            else
            {
                MPI_Bcast(&tempPivot,1,MPI_FLOAT,processId,Current_Comm); //SECOND BROADCAST : SENDING PIVOT   k ton stelnw sto lao
                pivot=tempPivot;
            }
        }
        
        if(stillActive==1)   //an eksakolouthw na eimai active, trexw tin partition.. k to count kommati to opio eimape kapou exei problima
        {
            partition(arrayToUse,elements,pivot,&arraySmall,&arrayBig,&endSmall,&endBig);
        }
        else
        {
            endBig=0;
            endSmall=0;
        }
        //an eimai inactive stelnw endbig=0 gia to bigset pou den epireazei
        sumSets=0;
        MPI_Reduce(&endBig,&sumSets,1,MPI_INT,MPI_SUM,0,Current_Comm); //FIRST REDUCE : SUM OF BIG, stelnw ola ta bigset gia na athroistoun sotn master
        MPI_Bcast(&sumSets,1,MPI_INT,0,Current_Comm);
        if(oldSumSets==sumSets)
            checkIdentical=1;
        else
        {
            oldSumSets=sumSets;
            checkIdentical=0;
        }
        MPI_Bcast(&finalize,1,MPI_INT,0,Current_Comm);//an o master apo to keepbigset k apo to count apofasisei oti teleiwsame mou stelnei 1, alliws 0 sunexizoume
        if(finalize==1)
        {
            float median=0;
            validation(median,partLength,size,numberPart,processId,Current_Comm,countOuter);
            MPI_Barrier(Current_Comm);
            return ;
        }
        MPI_Bcast(&keepBigSet,1,MPI_INT,0,Current_Comm);//THIRD BROADCAST: Receive keepBigset boolean, edw lambanw an krataw to mikro i megalo set.
            //afou elaba ton keepbigset an eimai active krataw enan apo tous duo pinake small h big.. alliws den kanw tpt
            //edw antistoixa allazw tous pointers, k eksetazw an exw meinei xwris stoixeia tin opoia periptwsi sikwnw to dropoutflag k pio katw tha dilwsw na ginw inactive
        if(stillActive==1)
        {
            if(keepBigSet==1)
            {
                if(endBig==0)
                    dropoutflag=1;
                else
                {
                    arrayToUse=arrayBig;
                    elements=endBig;
                }
            }
            else if(keepBigSet==0)
            {
                if(endSmall==0)
                    dropoutflag=1;
                else
                {
                    arrayToUse=arraySmall;
                    elements=endSmall;
                }
            }
        }
        //edw einai ligo periploka grammeno, isws exei perita mesa alla, an eimai active k thelw na ginw inactive einai i prwti periptwsi, h deuteri einai eimai inactive hdh k i triti einai sunexizw dunamika
        if(dropoutflag==1 && stillActive==1)
        {
            MPI_Send(&dropoutflag,1,MPI_INT,0,1,Current_Comm); //FIRST SEND : send active or not;
            stillActive=0;
        }
        else if(stillActive==0)
        {
            dropoutflag=-1;
            MPI_Send(&dropoutflag,1,MPI_INT,0,1,Current_Comm); //FIRST SEND : send active or not;
        }
        else
        {
            dropoutflag=0;
            MPI_Send(&dropoutflag,1,MPI_INT,0,1,Current_Comm); //FIRST SEND : send active or not;
        }
    }
}

void transferPoints(float *distances,float median,floatType **pointsCoords,int partLength,int coordSize,
    int child_Id,int child_num,MPI_Comm Current_Comm,int *counts){
    int i,j,k,offset,*allDests,*dests,destSize,partnersCounter,*allCounters;
    if(child_Id==0)
        allCounters = (int*)malloc(child_num*sizeof(int));
    MPI_Status Stat;
    int myCounter = counts[0];
    int dest = (child_Id < child_num / 2) ? (child_Id + child_num / 2) : (child_Id - child_num / 2);
    #ifdef DEBUG
    printf("%d Sending counter to...%d, myCounter = [%d]\n",child_Id,dest,myCounter);
    #endif
    MPI_Gather(&myCounter, 1, MPI_INT, allCounters, 1, MPI_INT, 0,Current_Comm);
    //Process the counters array for proper scattering
    if(child_Id == 0){
        #ifdef DEBUG_TRANSFER
        printf("allCountersGathered = |");
        for(i = 0;i < child_num;i++)
            printf("%d|",allCounters[i]);
        printf("\n");
        #endif
        for(i = 0;i < child_num / 2;i++){
            j = allCounters[i];
            allCounters[i] = allCounters[i + child_num/2];
            allCounters[i + child_num/2] = j;
        }
        #ifdef DEBUG_TRANSFER
        printf("allCountersScattering = |");
        for(i = 0;i < child_num;i++)
            printf("%d|",allCounters[i]);
        printf("\n");
        #endif
    }
    MPI_Scatter(allCounters, 1, MPI_INT, &partnersCounter, 1, MPI_INT,0,Current_Comm);
    #ifdef DEBUG
        printf("%d Recieved counter from master, partnersCounter = [%d]\n",child_Id,partnersCounter);
    #endif

    MPI_Barrier(Current_Comm);
    for(i = 0,offset = 0;i < partLength && offset == 0;i++){
        if(((distances[i] > median) && (child_Id < dest)) || ((distances[i] <= median) && (child_Id > dest))){
            #ifdef DEBUG
                printf("%d Sending distances[%d] to...%d\n",child_Id,i,dest);
            #endif
            //pointBuffer = pointsCoords[i];
            MPI_Sendrecv_replace(pointsCoords[i], coordSize, MPI_floatType, dest, 1, dest, 1,Current_Comm, &Stat);
            //pointsCoords[i] = pointBuffer;
            //printf("%d Recieved from...%d, myCounter = [%d], partnersCounter = [%d],point = [%f,%f,%f]\n",child_Id,dest,myCounter,partnersCounter,pointsCoords[i][0],pointsCoords[i][1],pointsCoords[i][2]);
            if(child_Id < child_num / 2){
                myCounter--;
                partnersCounter++;
                if(myCounter == 0 || partnersCounter == partLength)
                    offset = 1;
            }
            else{
                myCounter++;
                partnersCounter--;
                if(myCounter == partLength || partnersCounter == 0)
                    offset = 1;
            }

        }
    }
    
    MPI_Gather(&myCounter, 1, MPI_INT, allCounters, 1, MPI_INT, 0,Current_Comm);
    MPI_Barrier(Current_Comm);
    if(child_Id==0){
        destSize = 0;
        #ifdef DEBUG_TRANSFER
            printf("Final allCounters = |");
        #endif
        for(j = 0;j < child_num;j++){
            if((j < child_num/2 && allCounters[j] > destSize) 
                || (j >= child_num/2 && partLength - allCounters[j] > destSize))
                destSize = (j < child_num/2) ? (allCounters[j]):(partLength - allCounters[j]);
            #ifdef DEBUG_TRANSFER
                printf("%d|",allCounters[j]);
            #endif
        }
        #ifdef DEBUG_TRANSFER
            printf("\n,destSize = %d\n",destSize);
        #endif
        allDests = (int*)malloc(destSize * child_num * sizeof(int));
        for(j = 0;j < destSize * child_num;j++)
            allDests[j] = -1;
        for(j = 0;j < child_num/2;j++){
            for(k = child_num/2;k < child_num && allCounters[j] > 0;k++){
                while(allCounters[k] < partLength && allCounters[j] > 0){
                    allCounters[k]++;
                    allCounters[j]--;
                    offset = 0; 
                    while(allDests[j*destSize + offset]!=-1){offset++;}
                    allDests[j*destSize + offset] = k;
                    offset = 0; 
                    while(allDests[k*destSize + offset]!=-1){offset++;}
                    allDests[k*destSize + offset] = j;
                }
            }            
        }
        #ifdef DEBUG_TRANSFER
        printf("allDests = ");
        for(j = 0;j < child_num;j++){
            printf("|");
            for(k = 0;k < destSize;k++)
                printf("%d,",allDests[j*destSize+k]);
        }
        printf("|\n");
        #endif
    }
    MPI_Bcast(&destSize,1,MPI_INT,0,Current_Comm);

    //Redistribute points that were not equally placed accross "symmetric" procceses
    if(destSize > 0){
        dests = (int*)malloc(destSize * sizeof(int));
        MPI_Barrier(Current_Comm);
        MPI_Scatter(allDests, destSize, MPI_INT, dests, destSize, MPI_INT,0,Current_Comm);
        //printf("I am %d my dests are %d\n",child_Id,dests[0]);    
        
        if(dests[0] != -1){
            for(offset = 0; //i's value has been left unchanged and the point transfer continues right from the point it has stopped
                i < partLength && offset < destSize && dests[offset]!=-1;
                i++){
                if(((distances[i] > median) && (child_Id < dest)) || ((distances[i] <= median) && (child_Id > dest))){
                    #ifdef DEBUG
                        printf("%d Sending distances[%d] to...%d\n",child_Id,i,dests[offset]);
                    #endif
                    //pointBuffer = pointsCoords[i];
                    MPI_Sendrecv_replace(pointsCoords[i], coordSize, MPI_floatType, dests[offset], 1, dests[offset], 1,Current_Comm, &Stat);
                    //pointsCoords[i] = pointBuffer;
                    offset++;
                }
                
            }

        }
        MPI_Barrier(Current_Comm);
    }
    #ifdef DEBUG_TRANSFER
        printf("I'm done %d, myCounter = %d, i = %d\n",child_Id,myCounter,i);
    #endif
    //return pointsCoords;
}

//Reoder an array according to a pivot value
void transferPointsST(float* distances,float median,floatType **pointsCoords,int size,int coordSize){
    int i,j;
    floatType *temp;
    //Initial value of j as size - 1 is implemented using the predecrement of j [--j]
    j = size-1;
    for(i = 0;i < size && i < j;i++){
        if(distances[i] > median){
            //Scan until you find a proper point to transfer
            while(distances[j] > median && i < j)
                j--;

            temp = pointsCoords[i];
            pointsCoords[i] = pointsCoords[j];
            pointsCoords[j] = temp;
            /*
            for(k = 0;k < coordSize;k++){
                temp = pointsCoords[i][k];
                pointsCoords[i][k] = pointsCoords[j][k];
                pointsCoords[j][k] = temp;
            }*/
            j--;
        }
    }
}

const char* getfield(char* line, int num)
{
    const char* tok;
    for (tok = strtok(line, ",");
            tok && *tok;
            tok = strtok(NULL, ",\n"))
    {
        if (!--num)
            return tok;
    }
    return NULL;
}

//Code used as in(using j < col in addition): https://gist.github.com/amirmasoudabdol/f1efda29760b97f16e0e
void read_csv(int row, int col, char *filename, floatType **data, int rowOffset, int colOffset){
        FILE *file;
        file = fopen(filename, "r");
        int j = 0;
        char line[4098];
        while ((j < rowOffset) && fgets(line, 4098, file))
            j++;

        int i = 0;
        while (fgets(line, 4098, file) && (i < row))
        {
            // double row[ssParams->nreal + 1];
            char* tmp = strdup(line);
    
            j = 0;
            const char* tok;
            for (tok = strtok(line, ","); tok && *tok && j < col + colOffset; j++, tok = strtok(NULL, ",\n"))
            {
                if(j >= colOffset){
                    data[i][j - colOffset] = atof(tok);
                    //printf("%f\t", data[i][j - colOffset]);
                }

                
            }
            //printf("\n");
            //printf("File Size = (%d,%d)\n",i,j);
            free(tmp);
            i++;
        }
        //printf("File Size = (%d,%d)\n",i,j);
    }
    void printPoints(floatType **pointsCoords,int size,int coordSize,int processId){
        int i, j;
        for(i=0; i<size; i++){
            printf("point %d = [ ",i+processId*size);
            for(j=0;j<coordSize;j++)
                printf("%f,",pointsCoords[i][j]);
        printf("\b]\n");
        }
    }

/*========================FIND MEDIAN FUNCTIONS====================================
 * ================================================================================
 * ================================================================================
*/


/****Partitions the Array into larger and smaller than the pivot values****/
void partition (float *array, int elements, float pivot, float **arraysmall, float **arraybig, int *endsmall, int *endbig)
{
    int right=elements-1;
    int left=0;
    int pos;
    if(elements==1)
    {
        if(pivot>array[0])
        {
            *endsmall=1;  //One value in the small part
            *endbig=0;   //Zero on the big one
            *arraysmall=array;   //There is no big array therefore NULL value
            *arraybig=NULL;
        }
        else if(pivot<=array[0])
        {
            *endsmall=0;    //The exact opposite of the above actions.
            *endbig=1;
            *arraysmall=NULL;
            *arraybig=array;
        }
    }
    else if(elements>1)
    {
        while(left<right)
        {
            while(array[left]<pivot)
            {
                left++;
                if(left>=elements)
                {
                    break;
                }
            }
            while(array[right]>=pivot)
            {
                right--;
                if(right<0)
                {
                    break;
                }
            }
            if(left<right)
            {
                swap_values(array,left,right);
            }
        }
        pos=right;
        if(pos<0)                   //Arrange the arrays so that they are split into two smaller ones.
        {                               //One containing the small ones. And one the big ones.
            *arraysmall=NULL;           //However these arrays are virtual meaning that we only save the pointer values of the beging and end
        }                               //of the "real" one.
        else
        {
            *arraysmall=array;
        }
        *endsmall=pos+1;
        *arraybig=&array[pos+1];
        *endbig=elements-pos-1;
    }
}


/***==============================================***/
/***==============================================***/
/***=============SERIAL SELECTION==============***/
/***==============================================***/
/***==============================================***/

float selection(float *array,int number)
{
    float *arraybig;
    float *arraysmall;
    int endsmall=0;
    int endbig=0;
    float *arraytobeused;
    int i;
    int counter=0;
    int k;
    float pivot;
    float median;
    k=(int)number/2+1;
    arraytobeused=array;
    for(;;)
    {
        pivot=arraytobeused[rand() % number];
        partition(arraytobeused,number,pivot,&arraysmall,&arraybig,&endsmall,&endbig);
        if(endbig>k)
        {
            number=endbig;
            arraytobeused=arraybig;
            for(i=0;i<endbig;i++)
            {
                if(pivot==arraybig[i])
                    counter++;
                else
                    break;
            }
            if(counter==endbig)
            {
                median=arraybig[0];
                break;
            }
            else
                counter=0;
            //end of count equals
        }
        else if(endbig<k)
        {
            number=endsmall;
            arraytobeused=arraysmall;
            k=k-endbig;
        }
        else
        {
            median=pivot;
            break;
        }
    }
    return median;
}

