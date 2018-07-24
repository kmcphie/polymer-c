#!/usr/bin/perl
# creates a number of job execution lines for submission to pbs
use strict;
use warnings;

my $seriesName = "CD3ZetaElectrostaticsPhosphorylation"; # command line
my $occupiediSitesFile = 'OccupiediSitesMouse.txt';
my $occupiediSitesFileNoSpace = 'OccupiediSitesMouseNoSpace.txt';

my $i0Max = 3;
my $i1Max = 3;

for (my $i0 = 1; $i0 <= $i0Max; $i0++)
{
	my $PARABOLADEPTH = 10**(0.5*($i0-1));

	for (my $i1 = 1; $i1 <= $i1Max; $i1++)
	{
		my $PARABOLAWIDTH = 10**(-1+0.5*($i1-1));
        
        
        open(my $fileToRead, '<', $occupiediSitesFile)
            or die "Could not open file '$occupiediSitesFile' $!";
        
        open(my $fileToReadNoSpace, '<', $occupiediSitesFileNoSpace)
            or die "Could not open file '$occupiediSitesFileNoSpace' $!";
        
        while (my $line = <$fileToRead>)
        {
            chomp $line;
            
            my $lineNoSpace = <$fileToReadNoSpace>;
            chomp $lineNoSpace;


            my $fileName = $seriesName . "ParabolaDepth" . "." . $i0 . "." . "ParabolaWidth" . "." . $i1 . "." . $. ;
            my $runName =  $fileName . ".pub";
	
            open (FOOD, ">pubs/$runName" );
            print FOOD << "EOF";

#!/bin/bash
#\$ -N $seriesName
#\$ -q free*,pub*,abio,bio
#\$ -ckpt blcr
#\$ -e logs/$runName.err
#\$ -o logs/$runName.log

cd /pub/laraclemens/polymer-c_runs/Mar212017ElectrostaticsPhosphorylation

echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`

# Run your executable and then stuff
./metropolis.out parameters.txt $fileName $line $lineNoSpace $PARABOLADEPTH $PARABOLAWIDTH

echo Finished at `date`

EOF
		close FOOD;
        }
        
        close $fileToRead;
        close $fileToReadNoSpace;
	}
}
