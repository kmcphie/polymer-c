/*** Allard Group jun.allard@uci.edu                    ***/

void binaryConversion();


/*******************************************************************************/
//  GLOBAL VARIABLES for binary conversion
/*******************************************************************************/

int leftover,base,num,no_of_1s;
double binary;

/*******************************************************************************/

void binaryConversion()
{
    
    /*
     * C program to accept a decimal number and convert it to binary
     * and count the number of 1's in the binary number
     */
    // Adapted from: http://www.sanfoundry.com/c-program-decimal-binary-count-1-binary/
    
    for (i=0;i<sizeOfRateMatrix;i++)
    {
        base     = 1;
        binary   = 0;
        no_of_1s = 0;
        num      = i;
        while (num > 0)
        {
            leftover = num % 2;
            /*  To count no.of 1s */
            if (leftover == 1)
            {
                no_of_1s++;
            }
            binary = binary + leftover * base;
            num = num / 2;
            base = base * 10;
        }
        binaryState[i] = binary;
        totalBound[i] = no_of_1s;
     }
    
    //debugging
    if(0)
    {
        // print matrix of number bound
        for (i=0;i<sizeOfRateMatrix;i++)
        {
            printf("Number bound: %d \n", totalBound[i]);
        }
        
        // print matrix of binary numbers
        for (i=0;i<sizeOfRateMatrix;i++)
        {
            printf("binaryState: %lf \n", binaryState[i]);
        }
    }
    
    
}



