#!/bin/bash

if [ -f stop ]
    then
    rm stop
fi


touch waiter
flag=1
for curr in $*
    do
		if [ $flag -eq 1 ]
            then
			flag=0
			if [ "$curr" == "found" ]
                then
                if [ -f commands ]
                    then
                    echo "removing previous 'commands'"
                    rm commands
                fi
                
		#		go.sh found/found
                #wait
                ritzab.sh found found states spam
                wait
                ritzas.sh found found states spam
                wait
	    		prevKrylovD.sh $curr $curr states
                wait
			fi
		else
			addStage.sh $prev $curr
            cp inc $curr/inc
            krylov1.sh $prev $curr states spam
            wait
            for f in `cat sums`
                do
                sum1.sh $curr $curr states $f
            done
            wait
            ritzb.sh $prev $curr states spam sums
            wait
            ritzs.sh $prev $curr states spam sums
            wait
         	prevKrylovD.sh  $curr $curr states
		fi
        prev=$curr
done
rm waiter
