# Parallelization Script
# July 13, 2016

NRequested=0            # initialize number of runs submitted

NRODS=113       # N=143 for CD3 in human and mouse

IRATIO=7        # Vary for different ligands (kinase, phosphotase, ZAP-70 etc)

BRATIO=5        # BRATIO = 5-7 for SH2 domain of ZAP-70

FORCE=0

VERBOSE=0

TESTRUN=4 # 0 = not test run, use first set of hardcoded iSites, 1 and 2 - use test run iSites

ISITELOCATION=-1

BSITELOCATION=-1

STIFFENRANGE=0 # -1 means don't stiffen

STIFFCASE=0 # 0 = not test run - CD3Zeta


ISITEFILE="iSites.txt"

BSITEFILE="bSites.txt"

BSITECOMMAND=2

WELLDEPTH=1 #literally no idea what this parameter should be - is it even positive?

DEBYE=0.3 #arbitrary number between 0 and 0.5 nm? so between 0 and 0.6 kuhn lengths?

RWALL=0 #arbitrary number greater than -15nm (size of bilayer) or -50 kuhn lengths? should RWALL be measured in kuhn lengths?

PHOSELECTRORANGE=0

ITERATIONS=1


#####################################################

TOTALITERATIONS=1 #for testing

#TOTALITERATIONS=`wc -l < OccupiediSitesMouse.txt`

echo "Length of file is $TOTALITERATIONS"

    # set number of runs submitted by checking the running processes for lines with the program name, set NRequested to number of lines (may include one more than actual number of runs, since it counts grep -c metropolis in its tally)

NRequested=`ps | grep -c metropolis`

# while number of iterations ran is less than or equal to total number of iterations desired, loop through runs


#mkdir CD3ZetaWellDepthSweepMembraneOnCatFiles


#for ((WELLDEPTH=0;WELLDEPTH<=50;WELLDEPTH=($WELLDEPTH+5)))
#do

#echo "WellDepth = $WELLDEPTH"

ITERATIONS=1

while (( $ITERATIONS <= $TOTALITERATIONS ))
    do

    # loop to periodically check how many runs are submitted

    while (( $NRequested >= $1 ))   # while number requested is greater than number of processors we want to use (user input in command line)

        do
            sleep 1     # wait 1 sec before checking again

            NRequested=`ps | grep -c metropolis` # check again

    done

        echo "Done sleeping."

    # loop to submit more runs until reach max number of processors we want to use ($1)

    while (( $NRequested < $1 && $ITERATIONS <= $TOTALITERATIONS ))
        do

###############Ignore these lines - used for Stiffening.  Currently still in input, so need to be read, but otherwise can be ignored.###############
            # read specified line of text file

            OCCUPIEDSITES="`awk 'NR==iter' iter=$ITERATIONS OccupiediSitesMouse.txt`"

            OCCUPIEDSITESNOSPACE="`awk 'NR==iter' iter=$ITERATIONS OccupiediSitesMouseNoSpace.txt`"

            # print to screen the line read
            echo "Line $ITERATIONS of file is $OCCUPIEDISITES"
################################

            # run program with specified parameters

            # ./metropolis.out ElectrostaticsMembraneWallTest $NRODS $IRATIO $BRATIO $FORCE $VERBOSE $TESTRUN $ISITELOCATION $BSITELOCATION $STIFFENRANGE $STIFFCASE "$OCCUPIEDSITES"  "$OCCUPIEDSITESNOSPACE" "$ISITEFILE" "$BSITEFILE" $BSITECOMMAND $WELLDEPTH $DEBYE $RWALL $PHOSELECTRORANGE &

            ./metropolis.out parameters.txt TestParameterInput -1 -1 $WELLDEPTH $DEBYE &

            # If user gives V or v as second command line argument, then code will be verbose. Any other input will result in non-verbose.
            if [[ $2 == "V" || $2 == "v" ]]
                then
                    # print to screen the process ID and the name of the run
                    echo "PID of MultipleBinding.$ITERATIONS is $!"
            fi

            # update number of running programs
            NRequested=`ps | grep -c metropolis`

            # increase iteration run by 1
            ITERATIONS=$(( ITERATIONS + 1 ))
    done
            echo "Done calling metropolis."
done

# wait for all background processes to finish before concatenating files
wait

echo "Done waiting for processes to finish."

# loop through all files, concatenate them into one file
#for ((IT=1; IT<=$TOTALITERATIONS; IT++))
#do

#cat CD3ZetaMembrane1WellDepth.$WELLDEPTH.$IT >> CD3ZetaMembrane1WellDepth.$WELLDEPTH.cat.txt

#done

#cp CD3ZetaMembrane1WellDepth.$WELLDEPTH.cat.txt ~/Documents/polymer-c_runs/Oct142016WellDepthSweepMembraneOn/CD3ZetaWellDepthSweepMembraneOnCatFiles

#mkdir CD3ZetaMembrane1WellDepth.$WELLDEPTH

#mv CD3ZetaMembrane1WellDepth.$WELLDEPTH.* ~/Documents/polymer-c_runs/Oct142016WellDepthSweepMembraneOn/CD3ZetaMembrane1WellDepth.$WELLDEPTH/

#done

