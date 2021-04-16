#!/bin/bash

if [ $# -eq 5 ]
    then
     src=$1
     dir=$2
     prev=$src/$krylovBuild.kry
     curr=$dir/$krylovSum.sum
     file=$dir/$ritzSolve.riz
    echo "*Body $file" > $file
    echo "*InputOutput" >> $file
    echo "read found/found" >> $file
    echo "read $dir/stage " >> $file
    for state in `cat $3`
        do
        echo "vector $prev.$state" >> $file
        for sum in `cat $5`
            do
            echo "vector $curr.$sum.$state" >> $file
        done
    done
    for h in `cat $4`
        do
        echo "sum $dir/$ritzBuild.$h.riz " >> $file
    done
    echo "control ritzPhase" >> $file
    echo ".InputOutput" >> $file
    echo "*Parameters" >>$file
    echo "  spam 0 " >> $file
    echo ".Parameters" >> $file
    echo ".Body" >> $file
    go.sh $file &
fi
