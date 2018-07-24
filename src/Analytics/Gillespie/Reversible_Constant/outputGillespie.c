/*** Allard Group jun.allard@uci.edu                    ***/

void outputGillespie();
void dataRecording();

/*******************************************************************************/
//  GLOBAL VARIABLES for output control
/*******************************************************************************/

double finalBinaryState,averagePhosphorylationFraction_Iter,averageBound_Iter, averageBound_End, averagePhosphorylationFraction_End;
double sum_phosphorylationFraction,sum_time,sum_bound, sum_time_End, sum_bound_End,sum_phosphorylationFraction_End;
int finalPhosphorylationState;
int endStorage_length;

/********************************************************************************************************/
void initialize_dataRecording()
{
    iter=0;
}

/********************************************************************************************************/
void outputGillespie()
{
    if (!verboseTF)
    {
        
        // initializations
        sum_bound = 0;
        sum_phosphorylationFraction = 0;
        sum_time = 0;
        

        // find average phosphorylation state over last 100 states
    //    for(i=0;i<numberStatesStored;i++)
    //    {
    //        sum_phosphorylationState += stateStorage[i];
    //    }
    //    averagePhosphorylationState = (double) sum_phosphorylationState/ (double) numberStatesStored;
    //
        
        
        // find average phosphorylation state over last 100 states, taking time spent in state into account
        for(i=0;i<numberStatesStored;i++)
        {
            sum_bound                   += (double) timeStorage[i] * (double) stateStorage[i];
            sum_phosphorylationFraction += (double) timeStorage[i] * ((double) stateStorage[i] / (double) iSiteTotal);
            
            sum_time                    += (double) timeStorage[i];
        }
        averageBound_Iter                    = (double) sum_bound / (double) sum_time;
        averagePhosphorylationFraction_Iter  = (double) sum_phosphorylationFraction/ (double) sum_time;
        
        
        
        // find average phosphorylation state over last x time, taking time spent in each state into account
        
        for(i=0;i<endStorage_length;i++)
        {
            sum_bound_End                   += (double) timeStorage_End[i] * (double) stateStorage_End[i];
            sum_phosphorylationFraction_End += (double) timeStorage_End[i] * ((double) stateStorage_End[i] / (double) iSiteTotal);
            
            sum_time_End                    += (double) timeStorage_End[i];
         
            if(0)
            {
                printf("timeStorage_End %d : %f \n", i, timeStorage_End[i]);
                printf("sum_time_End: %f \n", sum_time_End);
            }
            
        }
        
        if(sum_time_End<timeAvgDuration)
        {
            printf("Not Enough Data");
            exit(0);
        }
        averageBound_End                    = (double) sum_bound_End / (double) sum_time_End;
        averagePhosphorylationFraction_End  = (double) sum_phosphorylationFraction_End/ (double) sum_time_End;
        
        if(0)
        {
            printf("averageBound_End: %f \n", averageBound_End);
            printf("averagePhosphorylationFraction_End: %f \n", averagePhosphorylationFraction_End);
        }
        
        
        // look at final state
        finalPhosphorylationState = totalBound[finalState];
        finalBinaryState          = binaryState[finalState];
        
        /***************** OUTPUT ********************/

        outputFile = fopen(outputName, "a");

        
        fprintf(outputFile, "%e ", reverseRate);                            // 1
        fprintf(outputFile, "%d ", iSiteTotal);                             // 2
        fprintf(outputFile, "%f ", timeEnd);                                // 3
        fprintf(outputFile, "%f ", averagePhosphorylationFraction_End);     // 4
        fprintf(outputFile, "%f ", averageBound_End);                       // 5
        fprintf(outputFile, "%d ", endStorage_length);                      // 6
        fprintf(outputFile, "%f ", averagePhosphorylationFraction_Iter);    // 7
        fprintf(outputFile, "%f ", averageBound_Iter);                      // 8
        fprintf(outputFile, "%d ", finalState);                             // 9
        fprintf(outputFile, "%d ", finalPhosphorylationState);              // 10
        fprintf(outputFile, "%f ", finalBinaryState);                       // 11
        fprintf(outputFile, "%f ", finalTotalTime);                         // 12
        fprintf(outputFile, "%d ", it);                                     // 13
        fprintf(outputFile, "\n");

        fclose(outputFile);
    
    }

    
}


/********************************************************************************************************/

void dataRecording()
{
    
    
    /*********************************************************************/
    // Record last x time states
    if (timeTotal >= (timeEnd-timeAvgDuration))
    {
        timeStorage_End[iter]  = timeStep;
        stateStorage_End[iter] = totalBound[currentState];
        endStorage_length = iter+1;
        
        if(0)
        {
            fflush(stdout);
            printf("timeStorage_End iter %d : %f \n", iter, timeStorage_End[iter]);
            fflush(stdout);
            printf("stateStorage_End iter %d : %d \n", iter, stateStorage_End[iter]);
            
        }
        
        iter++;
        


    }

    /*********************************************************************/
    if (verboseTF)
    {
        
        if (it%100==0)
        {
            // initializations - clear each time
            sum_bound = 0;
            sum_phosphorylationFraction = 0;
            sum_time = 0;
            
            // find average phosphorylation state over last 100 states, taking time spent in state into account
            for(i=0;i<numberStatesStored;i++)
            {
                sum_bound                   += (double) timeStorage[i] * (double) stateStorage[i];
                sum_phosphorylationFraction += (double) timeStorage[i] * ((double) stateStorage[i] / (double) iSiteTotal);
                
                sum_time                    += (double) timeStorage[i];
            }
            averageBound_Iter                    = (double) sum_bound / (double) sum_time;
            averagePhosphorylationFraction_Iter  = (double) sum_phosphorylationFraction/ (double) sum_time;
        
            /***************** OUTPUT ********************/

            outputFile = fopen(outputName, "a");
            
            
            fprintf(outputFile, "%e ", reverseRate);                             // 1
            fprintf(outputFile, "%d ", iSiteTotal);                              // 2
            fprintf(outputFile, "%f ", timeEnd);                                 // 3
            fprintf(outputFile, "%f ", averagePhosphorylationFraction_Iter);     // 4
            fprintf(outputFile, "%f ", averageBound_Iter);                       // 5
            fprintf(outputFile, "%f ", timeTotal);                               // 6
            fprintf(outputFile, "%d ", it);                                      // 7
            fprintf(outputFile, "\n");
            
            fclose(outputFile);
        }

    }
    
    
    
    
    
    
    
    
}





