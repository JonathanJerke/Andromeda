#!/bin/bash

if [ -f stop ]
    then
    rm stop
fi

if [ -f summa ]
    then
    prev=$1
    curr=$2
    #dir prev
    #dir number
    addStage.sh $prev $curr
    cp $prev/inc $curr/inc
    cp summa $curr
    touch waiter
    ritzSB.sh  $1 $2
    wait
    ritzSC.sh  $2
    wait
    prevKrylovD.sh $2 $2 states
    wait

    rm waiter
    
    
else
    echo "make a list of commands in summa"
    exit
fi
