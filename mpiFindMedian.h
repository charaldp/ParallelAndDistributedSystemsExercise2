#define DEBUG
//#define DEBUG_DIST
//#define DEBUG_TRANSFER
#define DEBUG_MAIN
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

/***Calculate distance***/
void calculateDistances(float *distances,floatType **pointCoords,floatType *vantagePointCoords,int pointsLength,int cordSize);

/***Validates the stability of the operation****/
void validation(float median,int partLength,int size,float *numberPart,int processId,MPI_Comm Current_Comm);

/***Validates the stability of the operation (Single Threaded)****/
void validationST(float median,int size,float *numberPart,int processId);

/***Validates the stability of tranfer point transactions for parallel iteration***/
void validationPartition(float median,int size,float *numberPart,int processId,int noProcesses);

/***Validates the stability of tranfer point transactions for serial iteration***/
void validationPartitionST(float median,int size,float *numberPart,int l);

/****Part executed only by the Master Node****/
float masterPart(int noProcesses,int processId,int size,int partLength,float *numberPart,MPI_Comm Current_Comm); //MASTER NODE CODE

/***Executed only by Slave nodes!!*****/
void slavePart(int processId,int partLength,float *numberPart,int size,int noProcesses,MPI_Comm Current_Comm); //code here is for the cheap slaves :P

//Perform a transaction of points across an even number of processes
void transferPoints(float *distances,float median,floatType **pointsCoords,int partLength,int coordSize,
    int child_Id,int child_num,MPI_Comm Current_Comm,int count);

//Reoder an array according to a pivot value
void transferPointsST(float* distances,float* medians,floatType **pointsCoords,int size,int coordSize);

//Read floating point numbers from a csv file into a matrix <data>
void read_csv(int row, int col, char *filename, floatType **data, int rowOffset, int colOffset);

/****Partitions the Array into larger and smaller than the pivot values****/
void partition (float *array,int elements, float pivot, float **arraysmall, float **arraybig, int *endsmall, int *endbig);

/***Serial Selection***/
float selection(float *array,int number);

float* multiSelection(float *array,int size,int l);

#endif



