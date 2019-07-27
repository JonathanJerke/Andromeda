#!/bin/csh

if ( -f waiter || -f stop ) then
	foreach com ( $argv )
		#cat $com >>  $com.a
		if ( -f waiter ) then
		echo "$com" >> commands
		endif
	end	
else
	foreach com ( $argv )
		date > $com.hout
		cat $com | olympics >> $com.hout
		sleep 2
		date >> $com.hout
	end                                                                               
endif
