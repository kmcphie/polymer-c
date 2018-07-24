#!/usr/bin/perl
# creates a number of job execution lines for submission to pbs

my $seriesName = "CD3ZetaElectrostaticsRwall6"; # command line

my $i0Max = 13;
my $i1Max = 13;
my $i2Max = 64; # number of lines in file

for (my $i0 = 1; $i0 <= $i0Max; $i0++)
{
    for (my $i1 = 1; $i1 <= $i1Max; $i1++)
    {
        
        my $i2 = 1;
        while ( $i2 <= $i2Max )
        {

            my $runName = $seriesName . "WellDepth" . "." . $i0 . "Debye" . "." . $i1 . "." . $i2 . ".pub";
			system "qsub pubs/$runName";

			print $i2;
			$i2++;

        }
    }
}
