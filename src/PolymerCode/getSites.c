/*** Allard Group jun.allard@uci.edu                    ***/

void getSites();


/*******************************************************************************/
//  GLOBAL VARIABLES for output control
/*******************************************************************************/


/********************************************************************************************************/
void getSites()
{
    /********* INITIALIZE ISITES *******************/
    
   
    switch (iSiteInputMethod) //switch in Batch Script to specify which set of iSites you want to use
    {
            
        case 0:  // iSites initialized for human CD3Zeta-Chain

            for(nf=0;nf<NFil;nf++)
            {
                iSiteTotal[nf] = 6;

                for (iy=0;iy<iSiteTotal[nf];iy++) //initializes iSite array
                {
                    iSite[nf][iy]=0;
                }
            }
            
            for(nf=0;nf<NFil;nf++)
            {
                //specify iSite locations by rod number (located at n-1) (i.e. if you have a polymer of N=50, and want an iSite at the 36th rod, input iSite[]=35)

                iSite[nf][0]=20;
                iSite[nf][1]=31;
                iSite[nf][2]=59;
                iSite[nf][3]=71;
                iSite[nf][4]=90;
                iSite[nf][5]=101;
            }
            break;
            
        case 1: //do nothing, use command line input
            
            break;
            
        case 2: //input iSites from file
            
            iSiteList = fopen(iSiteFilename, "r");
            char line[200];
            iy=0;
            
            while (fgets(line, sizeof(line), iSiteList))
            {
                // Eventually allow different iSites for each filament
                for(nf=0;nf<NFil;nf++)
                {
                    iSite[nf][iy]=atoi(line);
                }
                // Eventually have different numbers of iSites for each filament
                iy++;
            }
            
            fclose(iSiteList);
            
            for(nf=0;nf<NFil;nf++)
            {
                iSiteTotal[nf]=iy;
            }

            break;
        
        case 3: // use last site as only iSite
        
            for(nf=0;nf<NFil;nf++)
            {
                iSiteTotal[nf] = 1;
                for(iy=0;iy<iSiteTotal[nf];iy++)
                {
                    iSite[nf][iy]=0;
                }
            }
            
            for(nf=0;nf<NFil;nf++)
            {
                iSite[nf][0] = N[nf]-1;
            }
        
            break;

    }
    
    //Warning for possible user error
    for(nf=0;nf<NFil;nf++)
    {
        for (iy=0; iy<iSiteTotal[nf];iy++)
        {
            if (iSite[nf][iy] >= N[nf])
            {
                printf("Warning! Site is located past end of polymer in filament %ld!",nf);
                fflush(stdout);
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
            
            for (iy=0;iy<iSiteTotal[nf];iy++)
            {
                printf("iSite: %ld\n", iSite[nf][iy]);
                fflush(stdout);
            }
            
            printf("iSiteTotal: %ld\n", iSiteTotal[nf]);
            fflush(stdout);
        }
    }
    
    /*****************************************************/
    /********* INITIALIZE BOUND ISITES *******************/
    /*****************************************************/
    
    if (MULTIPLE) //if looking at multiple binding (i.e. MULTIPLE set to 1 in driveM)
    {
        switch (bSiteInputMethod)
        {
            case 0:
                for(nf=0;nf<NFil;nf++)
                {
                    bSiteTotal[nf] = 4; //total number of iSites bound
                    
                    for (ib=0;ib<bSiteTotal[nf];ib++) //initializes iSite array
                    {
                        bSite[nf][ib]=0;
                    }
                }
                for(nf=0;nf<NFil;nf++)
                {
                    bSite[nf][0]=10;
                    bSite[nf][1]=10;
                    bSite[nf][2]=40;
                    bSite[nf][3]=40;
                }
                break;

            case 1: //do nothing, use command line input
                
                break;
                
            case 2: //bSites for multiple binding of ZAP-70 to CD3 Zeta mouse
                
                for(nf=0;nf<NFil;nf++)
                {
                    for (iy=0; iy<iSiteTotal[nf]; iy++)
                    {
                        iSiteOccupied[nf][iy]=0;
                    }
                }
                // eventually want to be able to have different occupancies for each filament
                for(nf=0;nf<NFil;nf++)
                {
                    sscanf(occupiedSites,"%lf_%lf_%lf_%lf_%lf_%lf", &iSiteOccupied[nf][0],&iSiteOccupied[nf][1],&iSiteOccupied[nf][2],&iSiteOccupied[nf][3], &iSiteOccupied[nf][4],&iSiteOccupied[nf][5]);
                }
                
                for(nf=0;nf<NFil;nf++)
                {
                    bSiteCounter=0;
                    for (iy=0;iy<iSiteTotal[nf];iy++)
                    {
                        if (iSiteOccupied[nf][iy]==1)
                        {
                            bSite[nf][bSiteCounter]=iSite[nf][iy];
                            bSiteCounter++;
                        }
                    }
                    bSiteTotal[nf]=bSiteCounter;
                }

                break;
                
            case 3: // read bound sites from file
                
                bSiteList = fopen(bSiteFilename, "r");
                char line[200];
                iy=0;
                
                while (fgets(line, sizeof(line), bSiteList))
                {
                    // eventually want to be able to have different bound ligands for each filament
                    for(nf=0;nf<NFil;nf++)
                    {
                        bSite[nf][iy]=atoi(line);
                    }
                    iy++;
                    
                }
                
                fclose(bSiteList);
                
                for(nf=0;nf<NFil;nf++)
                {
                    bSiteTotal[nf]=iy;
                }
                
                break;
                
            case 4: // use last site as only bSite
                
                for(nf=0;nf<NFil;nf++)
                {
                    bSiteTotal[nf] = 1;
                    for(ib=0;ib<bSiteTotal[nf];ib++)
                    {
                        bSite[nf][ib]=0;
                    }
                }
                for(nf=0;nf<NFil;nf++)
                {
                    bSite[nf][0] = N[nf]-1;
                }
                
                break;
                
        }
        
        //Warning for possible user error
        for(nf=0;nf<NFil;nf++)
        {
            for (iy=0;iy<bSiteTotal[nf];iy++)
            {
                if (bSite[nf][iy] >= N[nf])
                {
                    printf("Warning! Bound site is located past end of filament %ld!",nf);
                    fflush(stdout);
                }
            }
        }
        
        //for debugging - prints a list of the bSites
        if (TALKATIVE)
        {
            for(nf=0;nf<NFil;nf++)
            {
                printf("Filament: %ld \n", nf);
                fflush(stdout);
                
                for (ib=0; ib<bSiteTotal[nf]; ib++)
                {
                    printf("bSite[%ld][%ld] is %ld \n",nf, ib, bSite[nf][ib]);
                    fflush(stdout);
                }
                
                printf("bSiteTotal = %ld \n", bSiteTotal[nf]);
                fflush(stdout);
            }
        }

    }
    
    
    


}

/********************************************************************************************************/

