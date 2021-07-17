#!/bin/bash

if [ $# -eq 3 ]
    then
     src=$1
     dir=$2
     prev=$src/$ritzSolve.riz
    for i in `cat $3`
        do
        file=$dir/$krylovBuild.kry.$i
        echo "*Body $file" > $file
        echo "*InputOutput" >> $file
        echo "read found/found" >> $file
        echo "read $dir/stage " >> $file
        echo "vector $prev-$i" >> $file
        echo "control buildPhase ">> $file
        echo ".InputOutput" >> $file
        echo ".Body" >> $file
        go.sh $file
    done
    
    if [ -e include ] 
    then
    for i in `cat include`
        do
        file=$dir/$krylovBuild.kry.$i
        echo "*Body $file" > $file
        echo "*InputOutput" >> $file
        echo "read found/found" >> $file
        echo "read $dir/stage " >> $file
        echo "vector $prev-$i" >> $file
        echo "control buildPhase ">> $file
        echo ".InputOutput" >> $file
        echo ".Body" >> $file
        go.sh $file
    done
    fi
    
fi

