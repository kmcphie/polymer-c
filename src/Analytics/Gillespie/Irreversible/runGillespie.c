/*** Allard Group jun.allard@uci.edu                    ***/

void runGillespie();

void runGillespie()
{
    
    /****************************** Create State Matrix from File ***********************/

    sizeOfStateMatrix = (int) pow(2,iSiteTotal);
    
    
    char line[200];
    i=0;
    j=0;
    
    ratesFile = fopen(matrixName,"r");
    
    while (fgets(line, sizeof(line), ratesFile))
    {
        stateMatrix[i][j]=atof(line);
        j++;
        if ((j%iSiteTotal) == 0)
        {
            i++;
            j=0;
        }
            
    }
    
    fclose(ratesFile);
    
    //debugging
    if(0)
    {
        for (i=0;i<sizeOfStateMatrix;i++)
        {
            for (j=0;j<iSiteTotal;j++)
            {
                printf("%lf ", stateMatrix[i][j]);
            
            }
        
            printf("\n");
        }
    }

    
    /******************************* Gillespie ******************************************/
    
    
    timeSum=0;
    it=0;
    
    while (it < iterations && it < ITMAX)
    {
    
    
        timeTotal=0;
        
        currentState=0;
        stepCount = 0;
        path=0;
        
        while (currentState < pow(2,iSiteTotal)-1) //while loop vs for loop?
        {
            
            
            //initialize random time array and time step
            for (iy=0;iy<iSiteTotal;iy++)
            {
                randTime[iy]=0;
            }
            
            timeStep = INF;
            
            //Gillespie step
            
            for (iy=0;iy<iSiteTotal;iy++)
            {
                if (stateMatrix[currentState][iy]!=0)
                {
                    randTime[iy] = - log(TWISTER)/stateMatrix[currentState][iy]; //exponentially distributed random variable based on transition rate
                }
                else
                {
                    randTime[iy] = 0; //use 0 instead of infinity - then just remove these cases later
                }
            }
            
            //pick smallest of random times
            for (iy=0;iy<iSiteTotal;iy++)
            {
                if (randTime[iy]!= 0)  // 0 time is not an option
                {
                    if (randTime[iy]<timeStep)
                    {
                        timeStep = randTime[iy];
                        newState = iy;
                    }
                }
            }
            
            if (0)
            printf("This is the path chosen: %d\n", newState);
            
            // create path number
            //use this update for "forwards" transitionMatrix (i.e. forwards binary, backwards phosphorylation)
            path += (iSiteTotal-newState)*pow(10,(iSiteTotal-stepCount-1));

            
            // use this update for "backwards" transitionMatrix (i.e. backwards binary)
            //path += (newState+1)*pow(10,(iSiteTotal-stepCount-1));

            
            //update time and state
            timeTotal += timeStep;
            stepCount++;
            
            // use this update for "forwards" transitionMatrix (i.e. forwards binary, backwards phosphorylation)
            currentState += pow(2,newState);
            
            // use this update for "backwards" transitionMatrix
            //currentState += pow(2,iSiteTotal-newState-1);
            
            if (0)
            printf("Current State is: %d \n", currentState);
            
            
        }
        
        // record which path is used and how long it took
        pathArray[path][0]++;
        pathArray[path][1] += timeTotal;
        
        //for MFPT
        timeSum += timeTotal;
        
        
        if (TIME)
        {
            timeArray[it] = timeTotal;
        }
        
        it++;
    
    }
    
    outputGillespie();
    
    
}



