#!/bin/bash

if [ $# -eq 3 ]
    then
    src=$1
    dir=$2
    for h in `cat $3`
        do
            file=$dir/$ritzBuild.$h.riz
            echo "*Body $file " > $file
            echo "*InputOutput" >> $file
            echo "read found/found" >> $file
            echo "read $src/stage " >> $file
#            for s in `cat states`
 #               do
                    for i in `cat $dir/summa`
                        do
                        echo "vector $i" >> $file
                    done
  #          done
            echo "control ritzPhase" >> $file
            echo ".InputOutput" >> $file
            echo "*Parameters" >>$file
            echo "  spam $h " >> $file
            echo ".Parameters" >> $file
            echo ".Body" >> $file
            go.sh $file &
        done
fi
