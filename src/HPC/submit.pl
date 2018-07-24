#!/usr/bin/perl
# creates a number of job execution lines for submission to pbs

my $seriesName = "CD3ZetaElectrostatics"; # command line

my $i0Max = 5;
my $i1Max = 5;
my $i2Max = 5;
my $i3Max = 5;
my $i4Max = 5;


for (my $i0 = 1; $i0 <= $i0Max; $i0++)
{
    
    for (my $i1 = 1; $i1 <= $i1Max; $i1++)
    {
        
        for (my $i2 = 1; $i2 <= $i2Max; $i2++)
        {
            
            for (my $i3 = 1; $i3 <= $i3Max; $i3++)
            {
                
                my $i4 = 1;
                while ( $i4 <= $i4Max )
                {

                        my $runName = $seriesName . "PD" . "." . $i0 . "PW" . "." . $i1 . "WK" . "." . $i2 . "ER" . "." . $i3 . "ZR" . "." . $i4 . ".pub";
                        system "qsub pubs/$runName";

                        print $i4;
                        $i4++;

                }
            }
        }
    }
}
