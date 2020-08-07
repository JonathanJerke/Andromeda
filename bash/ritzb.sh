#!/bin/bash

if [ $# -eq 4 ]
    then
    src=$1
    dir=$2
    prev=$src/D.kry
    curr=$dir/b.sum
    for h in `cat $4`
        do
        file=$dir/C.$h.riz
        echo "*Body $file " > $file
        echo "*InputOutput" >> $file
        echo "read found/found" >> $file
        echo "read $dir/stage " >> $file
        for f  in `cat $3`
            do
            echo "vector $prev.$f" >> $file
            echo "vector $curr.$f" >> $file
        done
        echo "control ritzPhase" >> $file
        echo ".InputOutput" >> $file
        echo "*Parameters" >>$file
        echo "  spam $h " >> $file
        echo ".Parameters" >> $file
        echo ".Body" >> $file
        go.sh $file &
    done
fi
