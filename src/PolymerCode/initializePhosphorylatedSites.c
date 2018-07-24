/*** Allard Group jun.allard@uci.edu                    ***/

void initializePhosphorylatedSites();

/*******************************************************************************/
//  GLOBAL VARIABLES for output control
/*******************************************************************************/

int totalPhosphorylated[NFILMAX];

double phosiSites[NFILMAX][NMAX];
/*******************************************************************************/
void initializePhosphorylatedSites()
{
    //initializes phoshorylatedSites to 0 (none phosphorylated)
    for(nf=0;nf<NFil;nf++)
    {
        for(ty=0;ty<iSiteTotal[nf];ty++)
        {
            phosiSites[nf][ty]=0;
        }
    }
 
    if (TALKATIVE)
    {
        for(nf=0;nf<NFil;nf++)
        {
            printf("These are the occupied (phosphorylated) sites on filament %ld : %s\n", nf, occupiedSites);
        }
    }
    
    //read string and assign to double vector
    // 1 is occupied iSite (phosphorylated), 0 is unoccupied
    //read string or read and interpret which of CD3Zeta Mouse tyrosines are phosphorylated
    for(nf=0;nf<NFil;nf++)
    {
        // Eventually want to have different phosphorylated sites depending on filament
        sscanf(occupiedSites,"%lf_%lf_%lf_%lf_%lf_%lf", &phosiSites[nf][0],&phosiSites[nf][1],&phosiSites[nf][2],&phosiSites[nf][3], &phosiSites[nf][4],&phosiSites[nf][5]);
    }

    //print which of the iSites are phosphorylated
    if (TALKATIVE)
    {
        for(nf=0;nf<NFil;nf++)
        {
            for (iy=0;iy<iSiteTotal[nf];iy++)
            {
                printf("phosiSites[ %ld ][ %ld ] =  %f\n",nf,iy, phosiSites[nf][iy]);
                fflush(stdout);
            }
        }
    }

    //initializes phosphorylated residues to 0 (none phosphorylated)
    for(nf=0;nf<NFil;nf++)
    {
        for(i=0;i<N[nf];i++)
        {
            PhosphorylatedSites[nf][i] = 0;
        }
    }

    //Phosphorylate sites
    // Note this is ONLY correct if PhosElectroRange does not extend past end points of polymer
    for(nf=0;nf<NFil;nf++)
    {
        for(ty=0;ty<iSiteTotal[nf];ty++)
        {
            if(phosiSites[nf][ty]==1) //might want to check the truth value on this - equals for double?
            {
                // include exit in case phosphorylation goes too far (can fix this to work appropriately if needed in the future)
                if( (iSite[nf][ty]-PhosElectroRange < 0) || (iSite[nf][ty]+PhosElectroRange+1 > N[nf]) )
                {
                    printf("Error! Phosphorylation extends beyond polymer %ld !",nf);
                    fflush(stdout);
                    exit(0);
                }
                // assign which segments are phosphorylated
                for(i=(iSite[nf][ty]-PhosElectroRange);i<(iSite[nf][ty]+PhosElectroRange+1);i++)
                {
                    PhosphorylatedSites[nf][i]=1; //set that joint to "phosphorylated"
                }
            }
        }
    }

    // for debugging, count and print total segments "phosphorylated"
    if (TALKATIVE)
    {
        for(nf=0;nf<NFil;nf++)
        {
            totalPhosphorylated[nf] = 0;
            for (i=0;i<N[nf];i++)
            {
                if (PhosphorylatedSites[nf][i]==1)
                {
                    totalPhosphorylated[nf]++;
                }
            }
            printf("Total Phosphorylated on filament %ld : %d\n", nf,totalPhosphorylated[nf]);
            fflush(stdout);
        }
    }
}
/***********************************************************************************/

