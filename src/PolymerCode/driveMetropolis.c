/*** Allard Lab jun.allard@uci.edu                    ***/

#define TWISTER genrand_real3()
#define NFILMAX         10
#define NMAX            400
#define NTMAX           2e9
#define NTADAPT         20000
#define NTCHECK         200000
#define DCHIMIN         1e-4
#define NBINS           100
#define NBINSPOLYMER    3000
#define PI              3.14159265359
#define INF             1e14
#define DCHIINIT        0.1
#define KSCRITICAL      0.01
#define MEMBRANE        0
#define MULTIPLE        0
#define STIFFEN         0
#define ELECTRO         0
#define HARDWALL        0
#define BASEBOUND       0
#define CPMAX           1e8
#define TALKATIVE       1
#define VISUALIZE       1
#define CD3ZETA         0
#define BINDTRANSITION  0

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

/* General Global Variables */
char listName[100];
FILE *fList;
//
char paramsFilename[100], iSiteFilename[100], bSiteFilename[100], basicSiteFilename[100];
FILE *paramsFile, *iSiteList, *bSiteList, *basicSiteList;

long NFil,N[NFILMAX];
long iSite[NFILMAX][NMAX], iSiteTotal[NFILMAX], iSiteCurrent, iy,ty, stericOcclusion[NFILMAX][NMAX];
long Ncurrent;
double c0, c1, irLigand;

double ree[NFILMAX], rM[NFILMAX], rM2[NFILMAX], rMiSite[NFILMAX][NMAX], rM2iSite[NFILMAX][NMAX], rH[NFILMAX];

double reeFil[NFILMAX][NFILMAX];

long iseed;

double phi[NFILMAX][NMAX], theta[NFILMAX][NMAX], psi[NFILMAX][NMAX];
double phiPropose[NFILMAX][NMAX], thetaPropose[NFILMAX][NMAX], psiPropose[NFILMAX][NMAX];
double r[NFILMAX][NMAX][3],t[NFILMAX][NMAX][3], e1[NFILMAX][NMAX][3], e2[NFILMAX][NMAX][3],
       rBase[NFILMAX][3], tBase[NFILMAX][3], e1Base[NFILMAX][3], e2Base[NFILMAX][3],
       rPropose[NFILMAX][NMAX][3],tPropose[NFILMAX][NMAX][3], e1Propose[NFILMAX][NMAX][3], e2Propose[NFILMAX][NMAX][3];
double norm;
double iLigandCenter[NFILMAX][NMAX][3];

double RGlobal[3][3], RLocal[3][3];
double e1_dot_t, e2_dot_t, e2_dot_e1;

long st;
long nf, nf2, nfPropose;
long proposals[2], accepts[2], nt, iChi, i, iPropose, ix, iParam, ntNextStationarityCheck,i2, iStart;
double E, ENew, rate[2], dChi[2], dChiHere, Force;
long constraintProposalsTotal;

int iSiteInputMethod;
long commandiSites;
char *iSiteLocations;
char input[4*NMAX];
long j,m;
int verboseTF;

/* Convergence Global Variables */
double ksStatistic[NFILMAX];
long ntNextStationarityCheck, iBin;
int convergedTF, constraintSatisfiedTF;
long convergenceVariableCounts[NFILMAX][NBINS], convergenceVariableCountsPrevious[NFILMAX][NBINS];
long polymerLocationCounts[NFILMAX][NMAX][NBINSPOLYMER];

/* STIFFEN Global Variables */
double StiffenRange, StiffSites[NFILMAX][NMAX];
int stiffCase, totalStiff[NFILMAX];

char occupiedSites[4*NMAX],occupiedSitesNoSpace[NMAX];
double iSiteOccupied[NFILMAX][NMAX];

/* MULTIPLE Global Variables*/
int bSiteInputMethod;
double brLigand;
double bLigandCenter[NFILMAX][NMAX][3];
long bSite[NFILMAX][NMAX], bSiteTotal[NFILMAX], bSiteCurrent, ib, ib2;
long bSiteCounter;

double bLigandCenterPropose[NFILMAX][NMAX][3];

double deliveryDistance;
long stericOcclusionBase[NFILMAX];
long membraneOcclusion[NFILMAX][NMAX], membraneAndSegmentOcclusion[NFILMAX][NMAX];
double localConcCutoff;
int deliveryMethod;
long boundToBaseDeliver[NFILMAX][NMAX];

/* ELECTRO Global Variables */
double Eelectro, EelectroNew;
double Erepulsion, Zrepulsion;
double parabolaDepth, parabolaWidth, wallParabolaK;
double PhosphorylatedSites[NFILMAX][NMAX];
int PhosElectroRange;
long basicSite[NFILMAX][NMAX], BasicSitesYN[NFILMAX][NMAX], basicSiteTotal[NFILMAX], basicSiteCurrent, iBasic;

/* BASEBOUND Global Variables */
double baserLigand;
double baseLigandCenter[NFILMAX][3];
double baseCenter[3];

/* MULTIPLE FILAMENT Variables*/
double baseSepDistance;

/* Filament tail dimerization force */
double dimerForce;

/*******************************************************************************/
//  INCLUDES
/*******************************************************************************/

#include "outputControl.c"
#include "getParameters.c"
#include "getSites.c"
#include "initializeStiffSites.c"
#include "initializePhosphorylatedSites.c"
#include "getBasicSites.c"
#include "metropolisJoint.c"


/*******************************************************************************/
//  MAIN
/*******************************************************************************/

// arguments: listName,
int main( int argc, char *argv[] )
{

    if(argv[1]) //get name of parameters file
        strcpy(paramsFilename,argv[1]);
    if (TALKATIVE) printf("This is the parameter filename: %s\n", paramsFilename);

    // use text file to set parameters
    getParameters();

    // use rest of command line input to overwrite parameters as necessary
    if(argv[2])
    {
        if(atoi(argv[2])!=-1)
            strcpy(listName,argv[2]);
        if (TALKATIVE) printf("This is the output filename: %s\n", listName);
    }

    if(argv[3])
    {
        if(atoi(argv[3])!=-1)
            verboseTF = atoi(argv[3]);
        if (TALKATIVE) printf("This will print verbose: %d\n", verboseTF);
    }

    if(argv[4])
    {
        if(atoi(argv[4])!=-1)
            NFil = atoi(argv[4]);
        if (TALKATIVE) printf("This is the number of filaments: %ld\n", NFil);
    }

    if(argv[5])
    {
        if(atoi(argv[5])!=-1)
            baseSepDistance = atof(argv[5]);
        if (TALKATIVE) printf("This is the base separation distance: %lf\n", baseSepDistance);
    }

    if(argv[6])
    {
        if(atoi(argv[6])!=-1)
            Force = atof(argv[6]);
        if (TALKATIVE) printf("This is the force: %lf\n", Force);
    }
    
    if(argv[7])
    {
        if(atoi(argv[7])!=-1)
            dimerForce = atof(argv[7]);
        if (TALKATIVE) printf("This is the dimerization force: %lf\n", dimerForce);
    }
    
    
    /* Finish setting up initial variables */
    
    // assign Ntemp to each filament
    for(nf=0;nf<NFil;nf++)
    {
        N[nf]=Ntemp;
        if (TALKATIVE) printf("This is number of rods in filament %ld: %ld\n",nf, N[nf]);
    }

	iseed = RanInitReturnIseed(0);
	
	metropolisJoint();

	return 0;
	
} // finished main
