#!/bin/bash

if [ $# -eq 2 ]
    then
        src=$1
        file=$src.ana.$2
        echo "*Body $file" > $file
        echo "*InputOutput" >> $file
        echo "  read $src" >> $file
        echo ".InputOutput" >> $file
        echo "*Parameters ">> $file
        echo "  iterations 1">> $file
        echo "  canonRank $2">> $file
        echo ".Parameters" >> $file
        echo ".Body" >> $file
        go.sh $file
    else
        echo "vector-file ##canoncial-ranks"
fi
