//#define DEBUG
//#define DEBUG_DIST
//#define DEBUG_TRANSFER
//#define DEBUG_TRANSFER_ST
//#define DEBUG_MAIN
//#define DBBUG_VPST

#define floatType float
#define MPI_floatType MPI_FLOAT

MPI_Status Stat;
extern double vpTimeSum;
extern double knnTimeSum;

#ifndef MPIFINDMEDIAN_H_
#define MPIFINDMEDIAN_H_

/***Kills processes that have no values left in their arrays****/
void removeElement(int *array, int *size, int element);

/****Swaps two values in an array****/
void swap_values(float *array,int x,int y);

/*** Generate a random point set for each process ***/
void generatePoints(floatType **pointsCoords,int pointsLength,int coordSize,int cal);

/***Calculate distances of set of points from another point***/
void calculateDistances(float *distances,floatType **pointCoords,floatType *vantagePointCoords,int pointsLength,int cordSize);

/***Calculate Singular Distance between 2 points***/
float distance(floatType *point1Coords,floatType *point2Coords,int cordSize);

/***Calculate distances of subsets of points from a single point of each subset***/
void calculateDistancesST(float *distances,floatType **pointCoords,int *vantagePoints,
                            int pointsLength,int vantagePointsLength,int cordSize);
/***Validates the stability of the operation****/
void validation(float median,int partLength,int size,float *numberPart,int processId,MPI_Comm Current_Comm);

/***Validates the stability of the operation (Single Threaded)****/
void validationST(float *medians,int size,float *numberPart,int processId,int multiplicity);

/****Part executed only by the Master Node****/
float masterPart(int noProcesses,int processId,int size,int partLength,float *numberPart,MPI_Comm Current_Comm); //MASTER NODE CODE

/***Executed only by Slave nodes!!*****/
void slavePart(int processId,int partLength,float *numberPart,int size,int noProcesses,MPI_Comm Current_Comm); //code here is for the cheap slaves :P

//Perform a transaction of points across an even number of processes
void transferPoints(float *distances,float median,floatType **pointsCoords,int partLength,int coordSize,
    int child_Id,int child_num,MPI_Comm Current_Comm,int count);

//Reoder an array according to a pivot value
void transferPointsST(float* distances,float* medians,floatType **pointsCoords,int size,int coordSize);

/***Validates the stability of tranfer point transactions for parallel iteration***/
void validationPartition(float median,int size,float *numberPart,int processId,int noProcesses);

/***Validates the stability of tranfer point transactions for serial iteration***/
void validationPartitionST(float *medians,int size,float *numberPart,int l);

//Read floating point numbers from a csv file into a matrix <data>
void read_csv(int row, int col, char *filename, floatType **data, int rowOffset, int colOffset);

void printPoints(floatType **pointsCoords,int size,int coordSize,int processId);

/****Partitions the Array into larger and smaller than the pivot values****/
void partition (float *array,int elements, float pivot, float **arraysmall, float **arraybig, int *endsmall, int *endbig);

/***Serial Selection***/
float selection(float *array,int number);

float* multiSelection(float *array,int size,int l);

void knnValidation(floatType **pointsCoords,float *distances,int pointsLength,floatType *pointCoords,int k,int *globalIndicesFound,MPI_Comm Current_Comm);
#endif



