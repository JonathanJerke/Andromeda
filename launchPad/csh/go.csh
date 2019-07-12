#!/bin/csh

if ( -f waiter ) then
	foreach com ( $argv )
		#cat $com >>  $com.a
		echo "$com" >> commands
	end	
else 
	foreach com ( $argv )
		date > $com.hout
		cat $com | olympics >> $com.hout
		sleep 2
		date >> $com.hout
	end                                                                               
endif
