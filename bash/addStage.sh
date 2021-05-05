#!/bin/bash

if [ $# -eq 2 ]
    then
    dir=$2
    newStage=$dir/stage
    src=$1
    mkdir $dir
#NO NEED, FOR MODULARIZE    python3 $LAUNCH/../python3/touch.py $dir/T
    #if [ ! -f $newStage ]
     #    then
        echo "*Body print " > $newStage
        echo "*InputOutput" >> $newStage
        echo " read $src/stage ">> $newStage
        echo " read $dir/inc" >> $newStage
        cp inc $dir
        echo ".InputOutput" >> $newStage
        echo ".Body" >> $newStage
        cp $krylovMultiply $dir
        cp $krylovSum $dir
        cp $ritzBuild $dir
        cp $ritzSolve $dir
        cp $krylovBuild $dir
   # fi
fi

