
//#define DEBUG
//#define DEBUG_TRANSFER
//#define DEBUG_MAIN
#define floatType double
#define MPI_floatType MPI_DOUBLE

MPI_Status Stat;


#ifndef MPIFINDMEDIAN_H_
#define MPIFINDMEDIAN_H_

/***Kills processes that have no values left in their arrays****/
void removeElement(int *array, int *size, int element);

/****Swaps two values in an array****/
void swap_values(float *array,int x,int y);

/*** Generate a random point set for each process ***/
void generatePoints(floatType **pointsCoords,int pointsLength,int coordSize,int cal);

/****Calculate distance****/
void calculateDistances(float *distances,floatType **pointCoords,floatType *vantagePointCoords,int pointsLength,int cordSize);

/***Validates the stability of the operation****/
void validation(float median,int partLength,int size,float *numberPart,int processId,MPI_Comm Current_Comm,int *counts);

/***Validates the stability of the operation (Single Threaded)****/
void validationST(float median,int size,float *numberPart,int *counts,int processId);

/***Validates the stability of the Multi Thread Point Partition ***/
void validationPartition(float median,int size,float *numberPart,int *counts,int processId);

/****Part executed only by the Master Node****/
float masterPart(int noProcesses,int processId,int size,int partLength,float *numberPart,MPI_Comm Current_Comm,int *countOuter); //MASTER NODE CODE

/***Executed only by Slave nodes!!*****/
void slavePart(int processId,int partLength,float *numberPart,int size,int noProcesses,MPI_Comm Current_Comm,int *countOuter); //code here is for the cheap slaves :P

//Perform a transaction of points across an even number of processes
void transferPoints(float *distances,float median,floatType **pointsCoords,int partLength,int coordSize,
    int child_Id,int child_num,MPI_Comm Current_Comm,int *counts);

//Reoder an array according to a pivot value
void transferPointsST(float* distances,float median,floatType **pointsCoords,int size,int coordSize);


const char* getfield(char* line, int num);


void read_csv(int row, int col, char *filename, floatType **data, int rowOffset, int colOffset);

void partition (float *array,int elements, float pivot, float **arraysmall, float **arraybig, int *endsmall, int *endbig);

float selection(float *array,int number);

#endif /* _QSORT-SEQUENTIAL_H_ */



