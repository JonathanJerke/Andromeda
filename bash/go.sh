#!/bin/bash

if [ -f stop ]
    then
    exit
fi

if [ -f waiter ]
	then
		for com in $* 
			do
				echo $com >> commands
			done	
else
	for com in $*
		do	
            before=`date`
            before2=`date -d "$before" +%s `
            if [ -f noOverWrite ]
                then
                    if [ -f $com.hout ]
                        then
                            echo "passed over $com b/c noOverWrite"
                        else
                            date > $com.hout
                            cat $com | andromeda >> $com.hout
                            date >> $com.hout

                        fi
                else
                date > $com.hout
                cat $com | andromeda >> $com.hout
                date >> $com.hout
            fi
            after=`date`
            after2=`date -d "$after" +%s `
            seconds=`echo "$before2 $after2" | awk '{print $2 - $1}'`
			if [[ `grep FINIS $com.hout` =~ "FINIS." ]]
				then
                    echo "$com $seconds FINIS."
				else
					echo $com > stop
					echo "stopped on $com"
			fi	
		done
fi
