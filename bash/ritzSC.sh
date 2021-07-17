#!/bin/bash

if [ $# -eq 1 ]
    then
    dir=$1
     
    file=$dir/$ritzSolve.riz
    echo "*Body $file" > $file
    echo "*InputOutput" >> $file
    echo "read found/found" >> $file
    echo "read $dir/stage " >> $file
#    for s in `cat states`
#        do
            for i in `cat $dir/summa`
                do
                echo "vector $i" >> $file
            done
 #   done
    for h in `cat spam`
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
