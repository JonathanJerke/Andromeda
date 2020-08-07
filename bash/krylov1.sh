#!/bin/bash

if  [ $# -eq 4 ]
    then
    src=$1
    dir=$2
    prev=$src/D.kry
    for i in `cat $3`
        do
        for h in `cat $4`
            do
            file=$dir/B.kry.$h.$i
            echo "*Body $file" > $file
            echo "*InputOutput" >> $file
            echo "read found/found" >> $file
            echo "read $dir/stage " >> $file
            echo "vector $prev.$i" >> $file
            echo "control kry1Phase ">> $file
            echo ".InputOutput" >> $file
            echo "*Parameters" >> $file
            echo "spam $h" >> $file
            echo ".Parameters" >> $file
            echo ".Body" >> $file
            go.sh $file &
        done
    done
fi

