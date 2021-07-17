#!/bin/bash

if [ -f stop ]
    then
    rm stop
fi

    prev=$1
    curr=$2
    #dir prev
    #dir number
    addStage.sh $prev $curr
    cp inc $curr/inc
    touch waiter
    foundAgain.sh $curr
    wait         
    ritzbF.sh $prev $curr states spam sums
    wait
    ritzsF.sh $prev $curr states spam sums
    wait
    prevKrylovD.sh $2 $2 states
    wait

    rm waiter
    
