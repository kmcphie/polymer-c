/*** Allard Group jun.allard@uci.edu                    ***/

void runGillespie();
void initializeStoreStates();
void storeStates();

void runGillespie()
{
    
    /****************************** Create State Matrix from File ***********************/

    sizeOfRateMatrix = (int) pow(2,iSiteTotal);
    
    // import first rate matrix
    char line1[200];
    i=0;
    j=0;
    
    ratesFile = fopen(matrixName1,"r");
    
    while (fgets(line1, sizeof(line1), ratesFile))
    {
        rateMatrix1[i][j]=atof(line1);
        j++;
        if ((j%sizeOfRateMatrix) == 0)
        {
            i++;
            j=0;
        }
            
    }
    
    fclose(ratesFile);
    
    // import second rate matrix
    char line2[200];
    i=0;
    j=0;
    
    ratesFile = fopen(matrixName2,"r");
    
    while (fgets(line2, sizeof(line2), ratesFile))
    {
        // import reverse values, multiply by variable constant
        rateMatrix2[i][j]=reverseRate*atof(line2);
        j++;
        if ((j%sizeOfRateMatrix) == 0)
        {
            i++;
            j=0;
        }
        
    }
    
    fclose(ratesFile);
    
    // concatenate matrices into full rate matrix
    
    for (i=0;i<sizeOfRateMatrix;i++)
    {
        for (j=0;j<sizeOfRateMatrix;j++)
        {
            //efficient, but not as safe - include if/else/reject statements instead?
            rateMatrix[i][j] = rateMatrix1[i][j] + rateMatrix2[i][j];
        }
    }
            
    
    //debugging
    if(0)
    {
        
        
        // print rate matrix
        for (i=0;i<sizeOfRateMatrix;i++)
        {
            for (j=0;j<sizeOfRateMatrix;j++)
            {
                printf("%lf ", rateMatrix[i][j]);
            }
        
            printf("\n");
        }
    }
    
    /*************** Binary Conversion and Initialize Stored States ************************/
    
    binaryConversion();
    initializeStoreStates();
    
    /******************************* Gillespie ******************************************/
    
    
    it=0;
    timeTotal=0;
    currentState=0; //start at fully dephosphorylated
    
    initialize_dataRecording();
    
    // while less than number of desired steps or less than max steps
    while (timeTotal < timeEnd && it < ITMAX)
    {
        
        //initialize random time array and time step
        for (iy=0;iy<sizeOfRateMatrix;iy++)
        {
            randTime[iy]=0;
        }
        
        timeStep = INF;
        
        //Gillespie step
        
        for (iy=0;iy<sizeOfRateMatrix;iy++)
        {
            if (rateMatrix[currentState][iy]!=0)
            {
                randTime[iy] = - log(TWISTER)/rateMatrix[currentState][iy]; //exponentially distributed random variable based on transition rate
            }
            else
            {
                randTime[iy] = 0; //use 0 instead of infinity - then just remove these cases later
            }
        }
        
        //pick smallest of random times
        for (iy=0;iy<sizeOfRateMatrix;iy++)
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
    
        //debugging
        if (1)
        {
            if(abs(totalBound[currentState]-totalBound[newState])!=1)
            {
                printf("Illegal state change!");
                exit(0);
            }
        }
        
        //update time
        timeTotal += timeStep;
        
        // record current state and timestep before updating - matches state with time spent in that state
        /******************************* Store Last 100 States ******************************************/
        //update stored states
        storeStates();
        dataRecording();
        /*************************************************************************************************/
        
        
        
        //update state
        currentState = newState;
    
        
        if (0)
        printf("Current State is: %d \n", currentState);
        


        

        it++;
    
    }
    finalState = currentState;
    finalTotalTime = timeTotal;
    
    outputGillespie();
    
    
}

/*************************************Initialize Store States********************************************/
/************************************************************************************************************************/
void initializeStoreStates()
{
    //initialize state storage
    numberStatesStored = 100;
    
    for (i=0;i<numberStatesStored;i++)
    {
        stateStorage[i] = 0;
        timeStorage[i]  = 0;
    }
    
}
/******************************************Store States*************************************************************/
/************************************************************************************************************************/
void storeStates()
{
    for (i=0;i<numberStatesStored-1;i++)
    {
        stateStorage[i] = stateStorage[i+1];
        timeStorage[i]  = timeStorage[i+1];
    }
    
    stateStorage[numberStatesStored-1] = totalBound[currentState];
    timeStorage[numberStatesStored-1]  = timeStep;
    

    if(0)
    {
        for (i=0;i<numberStatesStored;i++)
        {
            printf("State Storage %d: %d\n", i, stateStorage[i]);
        }
    
        for (i=0;i<numberStatesStored;i++)
        {
            printf("Time Storage %d: %f\n", i, timeStorage[i]);
        }
    
    }
}
    


