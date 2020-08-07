#!/bin/bash

if [ $# -eq 2 ]
    then
    dir=$2
    newStage=$dir/stage
    src=$1
    mkdir $dir
    if [ ! -f $newStage ]
         then
        echo "*Body print " > $newStage
        echo "*InputOutput" >> $newStage
        echo " read $src/stage ">> $newStage
        echo " read $dir/inc" >> $newStage
        cp inc $dir
        echo ".InputOutput" >> $newStage
        echo ".Body" >> $newStage
    fi
fi

