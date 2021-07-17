#!/bin/bash

if  [ $# -eq 1 ]
    then
    dir=$1
    prev=$src/D.kry
           
            file=$dir/A.fnd
            echo "*Body $file" > $file
            echo "*InputOutput" >> $file
            echo "read found/found" >> $file
            echo "read $dir/stage " >> $file
            echo "control foundPhase ">> $file
            echo ".InputOutput" >> $file
            echo "*Parameters" >> $file
            echo "	foundation 100 " >> $file
	    echo "      levelFoundation 1 " >> $file
	    echo "      widthFoundation 0.5 " >> $file
            echo ".Parameters" >> $file
            echo ".Body" >> $file
            go.sh $file &

fi

