#!/bin/bash

if [ $# -eq 4 ]
    then
         src=$1
         dir=$2
         prev=$src/$krylovMultiply.kry
         for i in `cat $3`
            do
            file=$dir/$krylovSum.sum.$4.$i
            echo "*Body $file" > $file
            echo "*InputOutput" >> $file
            echo "read found/found" >> $file
            echo "read $dir/stage " >> $file
            for h in `cat $4`
                do
                echo "vector $prev.$h.$i" >> $file
            done
            echo "control sumPhase ">> $file
            echo ".InputOutput" >> $file
            echo ".Body" >> $file
            go.sh $file &
        done
fi
