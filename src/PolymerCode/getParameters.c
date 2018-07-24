/*** Allard Group jun.allard@uci.edu                    ***/

void getParameters();


/*******************************************************************************/
//  GLOBAL VARIABLES for output control
/*******************************************************************************/
char tmpString[100];
long Ntemp;

/********************************************************************************************************/
void getParameters()
{
    paramsFile = fopen(paramsFilename,"r");
    
    fscanf(paramsFile,"%s %s", tmpString, listName);
    if (TALKATIVE) printf("This is output file name: %s\n", listName);
    
    fscanf(paramsFile,"%s %ld", tmpString, &NFil);
    if (TALKATIVE) printf("This is number of filaments: %ld\n", NFil);
    
    // eventually want different N per filament
    fscanf(paramsFile,"%s %ld", tmpString, &Ntemp);
    if (TALKATIVE) printf("This is number of rods in each filament: %ld\n",Ntemp);

    fscanf(paramsFile,"%s %lf", tmpString, &baseSepDistance);
    if (TALKATIVE) printf("This is the filament base separation distance: %lf\n", baseSepDistance);
    
    fscanf(paramsFile,"%s %lf", tmpString, &irLigand);
    if (TALKATIVE) printf("This is ligand radius: %lf\n", irLigand);
    
    fscanf(paramsFile,"%s %lf", tmpString, &brLigand);
    if (TALKATIVE) printf("This is bound ligand radius: %lf\n", brLigand);
    
    fscanf(paramsFile,"%s %lf", tmpString, &baserLigand);
    if (TALKATIVE) printf("This is base bound ligand radius: %lf\n", baserLigand);
    
    fscanf(paramsFile,"%s %lf", tmpString, &Force);
    if (TALKATIVE) printf("This is force: %f\n", Force);
    
    fscanf(paramsFile,"%s %lf", tmpString, &dimerForce);
    if (TALKATIVE) printf("This is dimerization force: %f\n", dimerForce);
    
    fscanf(paramsFile,"%s %d", tmpString, &verboseTF);
    if (TALKATIVE) printf("This is verbose: %d\n", verboseTF);
    
    fscanf(paramsFile,"%s %d", tmpString, &iSiteInputMethod);
    if (TALKATIVE) printf("This is iSite Input Method: %d\n", iSiteInputMethod);
    
    // omit argument 8 & 9 from text file - easier to just read iSites in from text file for now
    
    fscanf(paramsFile,"%s %lf", tmpString, &StiffenRange);
    if (TALKATIVE) printf("This is stiffening range: %lf\n", StiffenRange);
    
    fscanf(paramsFile,"%s %d", tmpString, &stiffCase);
    if (TALKATIVE) printf("This is stiffness case: %d\n", stiffCase);
    
    fscanf(paramsFile,"%s %s", tmpString, occupiedSites);
    if (TALKATIVE) printf("This is occupied Sites: %s\n", occupiedSites);
    
    fscanf(paramsFile,"%s %s", tmpString, occupiedSitesNoSpace);
    if (TALKATIVE) printf("This is occupied Sites: %s\n", occupiedSitesNoSpace);
    
    fscanf(paramsFile,"%s %s", tmpString, iSiteFilename);
    if (TALKATIVE) printf("This is iSite filename: %s\n", iSiteFilename);
    
    fscanf(paramsFile,"%s %s", tmpString, bSiteFilename);
    if (TALKATIVE) printf("This is bSite filename: %s\n", bSiteFilename);
    
    fscanf(paramsFile,"%s %s", tmpString, basicSiteFilename);
    if (TALKATIVE) printf("This is basicSite filename: %s\n", basicSiteFilename);
    
    fscanf(paramsFile,"%s %d", tmpString, &bSiteInputMethod);
    if (TALKATIVE) printf("This is bSite input method: %d\n", bSiteInputMethod);
    
    fscanf(paramsFile,"%s %lf", tmpString, &parabolaDepth);
    if (TALKATIVE) printf("This is parabolaDepth: %lf\n", parabolaDepth);
    
    fscanf(paramsFile,"%s %lf", tmpString, &parabolaWidth);
    if (TALKATIVE) printf("This is parabolaCenter: %lf\n", parabolaWidth);
    
    fscanf(paramsFile,"%s %lf", tmpString, &wallParabolaK);
    if (TALKATIVE) printf("This is wall parabola K: %lf\n", wallParabolaK);
    
    fscanf(paramsFile,"%s %lf", tmpString, &Erepulsion);
    if (TALKATIVE) printf("This is Erepulsion: %lf\n", Erepulsion);
    
    fscanf(paramsFile,"%s %lf", tmpString, &Zrepulsion);
    if (TALKATIVE) printf("This is Zrepulsion: %lf\n", Zrepulsion);
    
    fscanf(paramsFile,"%s %d", tmpString, &PhosElectroRange);
    if (TALKATIVE) printf("This is PhosElectroRange: %d\n", PhosElectroRange);
    
    fscanf(paramsFile,"%s %lf", tmpString, &localConcCutoff);
    if (TALKATIVE) printf("This is localConcCutoff: %lf\n", localConcCutoff);
    
    fclose(paramsFile);
    
}

/********************************************************************************************************/

