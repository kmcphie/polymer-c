/*** Allard Group jun.allard@uci.edu                    ***/

void getBasicSites();


/*******************************************************************************/
//  GLOBAL VARIABLES for output control
/*******************************************************************************/


/********************************************************************************************************/
void getBasicSites()
{
    /********* INITIALIZE BASIC AND TYROSINE SITES *******************/
    

    basicSiteList = fopen(basicSiteFilename, "r");
    char line[200];
    iBasic=0;
    
    while (fgets(line, sizeof(line), basicSiteList))
    {
        // For now, let all filaments have the same basic sites
        // Eventually allow each filament to have different basic sites
        for(nf=0;nf<NFil;nf++)
        {
            basicSite[nf][iBasic]=atoi(line);
        }
        // Eventually count sites separately for each filament
        iBasic++; // count how many basic sites per filament
    }
    
    fclose(iSiteList);
    
    for(nf=0;nf<NFil;nf++)
    {
        basicSiteTotal[nf]=iBasic;
    }

    for(nf=0;nf<NFil;nf++)
    {
        for (iBasic=0; iBasic<basicSiteTotal[nf];iBasic++)
        {
            if (basicSite[nf][iBasic] >= N[nf])
            {
                printf("Warning! Site is located past end of polymer number %ld !",nf);
                fflush(stdout);
            }
        }
    }
    
    //initializes basic residues to 0
    for(nf=0;nf<NFil;nf++)
    {
        for(i=0;i<N[nf];i++)
        {
            BasicSitesYN[nf][i] = 0;
        }
    }
    
    //Full residue list of basic yes or no
    for(nf=0;nf<NFil;nf++)
    {
        for (i=0;i<N[nf];i++)
        {
            for(iBasic=0;iBasic<basicSiteTotal[nf];iBasic++)
            {
                if(i==basicSite[nf][iBasic])
                {

                    BasicSitesYN[nf][i]=1; //set that joint to "basic"

                }
            }
        }
    }
    
    //for debugging - prints a list of the iSites
    
    if (TALKATIVE)
    {
        for(nf=0;nf<NFil;nf++)
        {
            printf("Filament: %ld\n", nf);
            fflush(stdout);
            
            for (iBasic=0;iBasic<basicSiteTotal[nf];iBasic++)
            {
                printf("basicSite: %ld\n", basicSite[nf][iBasic]);
                fflush(stdout);
            }
            
            printf("basicSiteTotal: %ld\n", basicSiteTotal[nf]);
            fflush(stdout);
        }
    }


}

/********************************************************************************************************/

