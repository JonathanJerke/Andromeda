#!/bin/bash

export LAUNCHER_JOB_FILE=launcher
echo "`date`"

if [ -f hosts ]
    then
        export LAUNCHER_NHOSTS=`cat hosts`
    else
        echo "set hosts for speedup"
        export LAUNCHER_NHOSTS=1
    fi

export LAUNCHER_WORKDIR=`pwd`

COM1=`head -n 1 launcher.com`
HEADER=`header.sh $COM1`
if [ -f $HEADER ]
    then
        export LAUNCHER_PPN=`cat $HEADER`
    elif [ -f ppn ]
        then
            export LAUNCHER_PPN=`cat ppn`
    else
        echo "set ppn for speedup"
        export LAUNCHER_PPN=4
fi

if [ -f jobDone ]
    then
        rm jobDone
    fi

echo "touch jobDone" >> launcher

for i in `seq 1 100`
	do
		$LAUNCHER_DIR/paramrun
		if [ -f jobDone ]
            then
                exit
            else
                echo "retry to connect with TACC-launcher"
                sleep 100
		fi
	done
