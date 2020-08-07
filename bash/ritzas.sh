#!/bin/bash

if [ $# -eq 4 ]
    then
     src=$1
     dir=$2
     prev=$src/D.kry
     curr=$dir/b.sum
     file=$dir/c.riz
    echo "*Body $file" > $file
    echo "*InputOutput" >> $file
    echo "read found/found" >> $file
    echo "read $dir/stage " >> $file
    echo "vector found/found" >> $file
    for h in `cat $4`
        do
        echo "sum $dir/C.$h.riz " >> $file
    done
    echo "control ritzPhase" >> $file
    echo ".InputOutput" >> $file
    echo "*Parameters" >>$file
    echo "  spam 0 " >> $file
    echo ".Parameters" >> $file
    echo ".Body" >> $file
    go.sh $file &
fi
