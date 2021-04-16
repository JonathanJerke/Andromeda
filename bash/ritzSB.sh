#!/bin/bash

if [ $# -eq 2 ]
    then
    src=$1
    dir=$2
    for h in `cat spam`
        do
            file=$dir/$ritzBuild.$h.riz
            echo "*Body $file " > $file
            echo "*InputOutput" >> $file
            echo "read found/found" >> $file
            echo "read $src/stage " >> $file
            for s in `cat states`
                do
                    for i in `grep $krylovBuild.kry.$s $dir/summa`
                        do
                        echo "vector $i" >> $file
                    done
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
