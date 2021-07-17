#!/bin/bash

if [ -f stop ] 
	then
        echo "rm'd stop"
        rm stop
	fi

if [ -f waiter ] 
	then
        echo "rm'd waiter"
        rm waiter
	fi

if [ -f launcher ] 
	then
        rm launcher
        rm launcher.com
	fi

if [ $# -eq 0 ]
	then
	com0=`head -n 1 commands`
	h=`header.sh $com0`
fi

if [ $# -gt 0 ]
	then
	h=$1
    echo "set to start at $h"
fi

if [ $# -gt 1 ]
    then
    f=$2
    echo "set to end before $f"
else
    f="x"
fi

if [ -f spn ]
then
    export LAUNCHER_SPN=`cat spn`
else
    echo "set Single Process Number for speedup"
    export LAUNCHER_SPN=100000000
fi


for com in `cat commands`
    do
        if [ -f stop ]
        then
            exit
        fi

        ch=`header.sh $com`
        if [ -f launcher ]
            then
                wc=`wc -l launcher`
                wc=`header.sh $wc`
            else
                wc=0
        fi
        if [ "$h" == "$ch" ]
            then
                echo "go.sh $com " >> launcher
                echo $com >> launcher.com

        elif [ $wc -gt 0 ]
            then
            if [ $wc -gt $LAUNCHER_SPN ]
                then
                before=`date`
                before2=`date -d "$before" +%s `
                htc.sh
                after=`date`
                after2=`date -d "$after" +%s `
                seconds=`echo "$before2 $after2" | awk '{print $2 - $1}'`
                echo "summary htc $h $seconds FINIS."
            else
                before=`date`
                before2=`date -d "$before" +%s `
                go.sh `cat launcher.com`
                after=`date`
                after2=`date -d "$after" +%s `
                seconds=`echo "$before2 $after2" | awk '{print $2 - $1}'`
                echo "summary go $h $seconds FINIS."
            fi
            rm launcher
            rm launcher.com
            
            h=$ch
            if [[ "$h" =~ "$f" ]]
                then
                    echo "FINIS."
                    exit
                fi
                
                
            echo "go.sh $com " >> launcher
            echo $com >> launcher.com
        fi
done
    
if [ -f stop ]
    then
        exit
    fi

    
if [ -f launcher ]
    then
        wc=`wc -l launcher`
        wc=`header.sh $wc`
    else
        wc=0
fi
if [ $wc -gt 0 ]
    then
    if [ $wc -gt $LAUNCHER_SPN ]
    then
        before=`date`
        before2=`date -d "$before" +%s `
        htc.sh
        after=`date`
        after2=`date -d "$after" +%s `
        seconds=`echo "$before2 $after2" | awk '{print $2 - $1}'`
        echo "summary htc $h $seconds FINIS."
    else
        before=`date`
        before2=`date -d "$before" +%s `
        go.sh `cat launcher.com`
        after=`date`
        after2=`date -d "$after" +%s `
        seconds=`echo "$before2 $after2" | awk '{print $2 - $1}'`
        echo "summary go $h $seconds FINIS."
    fi
else
    echo "null parsed commands 'launcher'"
    exit
fi
