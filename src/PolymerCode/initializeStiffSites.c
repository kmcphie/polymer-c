/*** Allard Group jun.allard@uci.edu                    ***/

void initializeStiffSites();

/*******************************************************************************/
//  GLOBAL VARIABLES for initializing stiff sites
/*******************************************************************************/
int stiffEnd, stiffStart;
double stiffiSites[NFILMAX][NMAX];
//
/********************************************************************************************************/
void initializeStiffSites()
{
    //initializes stiffiSites to 0 (none phosphorylated)
    for(nf=0;nf<NFil;nf++)
    {
        for(ty=0;ty<iSiteTotal[nf];ty++)
        {
            stiffiSites[nf][ty]=0;
        }
    }

    if (TALKATIVE)
    {
        for(nf=0;nf<NFil;nf++)
        {
            printf("These are the occupied sites for filament %ld: %s\n", nf, occupiedSites);
        }
    }
    
    //read string and assign to double vector
    // 1 is occupied iSite (phosphorylated), 0 is unoccupied
    switch (stiffCase)
    {
        case 0:
            
            for(nf=0;nf<NFil;nf++)
            {
                // eventually want to be able have different phosphorylation settings per iSite per filament
                sscanf(occupiedSites,"%lf_%lf_%lf_%lf_%lf_%lf", &stiffiSites[nf][0],&stiffiSites[nf][1],&stiffiSites[nf][2],&stiffiSites[nf][3],&stiffiSites[nf][4],&stiffiSites[nf][5]);
            }
            break;
            
            // could include a case to read occupiedSites from a file of either locations or of 0,1s
            // for now, this is only useful for CD3Zeta parameters (anything with six iSites)
            // could get rid of switch altogether for now
    }
    
    // for debugging, print which iSites are declared stiff
    if (TALKATIVE)
    {
        for(nf=0;nf<NFil;nf++)
        {
            for (iy=0;iy<iSiteTotal[nf];iy++)
            {
                printf("stiffiSites[ %ld ][ %ld ] =  %f\n",nf,iy, stiffiSites[nf][iy]);
                fflush(stdout);
            }
        }
    }
    
    //initializes stiffened rods to 0 (none stiff)
    for(nf=0;nf<NFil;nf++)
    {
        for(i=0;i<N[nf];i++)
        {
            StiffSites[nf][i] =0;
        }
    }

    /********************************************************/
    /******************* STIFFEN SEGMENTS *******************/
    /********************************************************/
    
    for(nf=0;nf<NFil;nf++)
    {
        for(ty=0;ty<iSiteTotal[nf];ty++)
        {
            if(stiffiSites[nf][ty]==1) //might want to check the truth value on this - equals for double?
            {
                // set beginning of stiffening range
                if(iSite[nf][ty]-StiffenRange >= 0)
                {
                    stiffStart=iSite[nf][ty]-StiffenRange;
                }
                else // if stiffenrange goes below 0, start stiffening at 0
                {
                    stiffStart=0;
                }
                
                //set end of stiffening range
                // if stiffenrange goes above N, end at N
                if(iSite[nf][ty]+StiffenRange+1 >= N[nf])
                {
                    stiffEnd=N[nf];
                }
                else
                {
                    stiffEnd=iSite[nf][ty]+StiffenRange+1;
                }
                
                // declare which segments are stiff, exclusive of right endpoint (because included +1 above)
                for(i=stiffStart;i<stiffEnd;i++)
                {
                    StiffSites[nf][i]=1; //set that joint to "stiff"
                }
            }
        }
    }

    // for debugging, count and print the total number of stiff segments
    if (TALKATIVE)
    {
        for(nf=0;nf<NFil;nf++)
        {
            totalStiff[nf] = 0;
            for (i=0;i<N[nf];i++)
            {
                if (StiffSites[nf][i]==1)
                {
                    totalStiff[nf]++;
                }
            }
            printf("Total Stiff on filament %ld: %d\n",nf, totalStiff[nf]);
            fflush(stdout);
            
            if (totalStiff[nf] >= N[nf])
            {
                // Include error for completely stiff
                // May cause convergence problems?
                printf("Error! Filament %ld is completely stiff!\n",nf);
                fflush(stdout);
                
                exit(0);
            }
        }
    }
}

/********************************************************************************************************/

