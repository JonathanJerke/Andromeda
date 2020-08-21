#!/bin/bash

if [ $# -eq 2 ]
    then
    mkdir $1
    cd $1
    mkdir found
    cp $LAUNCH/control/found found
    cp $LAUNCH/control/stage found
    cp $LAUNCH/control/two .
    cp $LAUNCH/control/inc .
    cp $LAUNCH/symmetry/$2 symmetry
    echo "spam" > sums
##
if [ -f ../cnfg ]
    then
    line=`grep states ../cnfg`
    states=`header2.sh $line`
    seq 1 $states  > states
fi

if [ -f ../cnfg ]
    then
    line=`grep spam ../cnfg`
    spam=`header2.sh $line`
    seq 1 $spam  > spam
fi


else
    echo "start.sh <name> SA"
fi

