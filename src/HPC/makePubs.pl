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
	my $PARABOLADEPTH = 10**(-2+($i0-1));

	for (my $i1 = 1; $i1 <= $i1Max; $i1++)
	{
		my $PARABOLAWIDTH = 10**(-2+($i1-1));
        
        for (my $i2 = 1; $i2 <= $i2Max; $i2++)
        {
            my $WALLPARABOLAK = 10**(-2+($i2-1));
            
            for (my $i3 = 1; $i3 <= $i3Max; $i3++)
            {
                my $EREPULSION = 10**(-2+($i3-1));
                
                for (my $i4 = 1; $i4 <= $i4Max; $i4++)
                {
                    my $ZREPULSION = 10**(-2+($i4-1));

        my $fileName = $seriesName . "PD" . "." . $i0 . "PW" . "." . $i1 . "WK" . "." . $i2 . "ER" . "." . $i3 . "ZR" . "." . $i4;
		my $runName =  $fileName . ".pub";
	
		open (FOOD, ">pubs/$runName" );
		print FOOD << "EOF";

#!/bin/bash
#\$ -N $seriesName
#\$ -q free*,pub*,abio,bio
#\$ -ckpt blcr
#\$ -e logs/$runName.err
#\$ -o logs/$runName.log

cd /pub/lclemens/polymer-c_runs/20171121CD3ZetaElectrostaticsSweep1/Dephos

echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`

# Run your executable and then stuff
./metropolis.out parameters.txt $fileName -1 -1 $PARABOLADEPTH $PARABOLAWIDTH $WALLPARABOLAK $EREPULSION $ZREPULSION

echo Finished at `date`

EOF
		close FOOD;	
	}	
}
