/*** Allard Lab lclemens@uci.edu                   ***/

#define TWISTER genrand_real3()
#define NMAX       400
#define NTMAX      1e9
#define NTADAPT    20000
#define NTCHECK    200000
#define DCHIMIN    1e-4
#define NBINS      100
#define PI         3.14159265359
#define INF        1e14
#define ISITEMAX   9
#define STATEMAX   1000
#define ITMAX      1000000000
#define TALKATIVE  1

#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include "twister.c"


/*******************************************************************************/
//  GLOBAL VARIABLES
/*******************************************************************************/

char matrixName1[1000];
char matrixName2[1000];
FILE *ratesFile;
long iseed;

char outputName[1000];
FILE *outputFile;


double timeTotal,randTime[STATEMAX],timeStep,timeEnd;
int currentState,iy,it,iterations;

double rateMatrix[STATEMAX][STATEMAX],rateMatrix1[STATEMAX][STATEMAX],rateMatrix2[STATEMAX][STATEMAX];
int i,j,k,iter;
int iSiteTotal,newState;

int sizeOfRateMatrix;
int totalBound[STATEMAX];
double binaryState[STATEMAX];
int verbose, summaryOn;
int stateStorage[100000],numberStatesStored;
double timeStorage[100000];

double timeAvgDuration;
int stateStorage_End[ITMAX];
double timeStorage_End[ITMAX];

double reverseRate;

double finalTotalTime;
int finalState;
int verboseTF;

/*******************************************************************************/
//  INCLUDES
/*******************************************************************************/

#include "binaryConversion.c"
#include "outputGillespie.c"
#include "runGillespie.c"


/*******************************************************************************/
//  MAIN
/*******************************************************************************/

// arguments: listName, dimension, nBar, N, dt
int main( int argc, char *argv[] )
{

    printf("This program is starting.\n");
    
    if(argv[1]) // matrixName
        strcpy(matrixName1, argv[1]);
    if (TALKATIVE) printf("This is argument 1: %s\n", matrixName1);
    
    
    if(argv[2]) // matrixName
        strcpy(matrixName2, argv[2]);
    if (TALKATIVE) printf("This is argument 2: %s\n", matrixName2);
    
    if(argv[3]) //iSiteTotal
        iSiteTotal = atoi(argv[3]);
    if (TALKATIVE) printf("This is argument 3: %d\n", iSiteTotal);
    
    if(argv[4]) //total number of iterations
        timeEnd = atof(argv[4]);
    if (TALKATIVE) printf("This is argument 4: %f\n", timeEnd);
    
    if(argv[5]) //total number of iterations
        reverseRate = atof(argv[5]);
    if (TALKATIVE) printf("This is argument 5: %f\n", reverseRate);
    
    if(argv[6]) //output file name
        strcpy(outputName, argv[6]);
    if (TALKATIVE) printf("This is argument 6: %s\n", outputName);
    
    if(argv[7]) //verbose T/F
        verboseTF = atoi(argv[7]);
    if (TALKATIVE) printf("This is argument 7: %d\n", verboseTF);
    
    if(argv[8]) //time to average over
        timeAvgDuration = atof(argv[8]);
    if (TALKATIVE) printf("This is argument 8: %f\n", timeAvgDuration);
    
    
	iseed = RanInitReturnIseed(0);
	
	runGillespie();

	return 0;
	
} // finished main
