#!/bin/bash

if [ $# -eq 3 ]
    then
     src=$1
     dir=$2
     prev=$src/c.riz
    for i in `cat $3`
        do
        file=$dir/D.kry.${i}
        echo "*Body $file" > $file
        echo "*InputOutput" >> $file
        echo "read found/found" >> $file
        echo "read $dir/stage " >> $file
        echo "vector $prev-${i}" >> $file
        echo "control buildPhase ">> $file
        echo ".InputOutput" >> $file
        echo ".Body" >> $file
        go.sh $file
    done
fi

