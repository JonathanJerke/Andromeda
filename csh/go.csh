#!/bin/csh

if ( -f waiter || -f stop ) then
	foreach com ( $argv )
		if ( -f waiter ) then
			echo "$com" >> commands
		endif
	end	
else

 
	foreach com ( $argv )	
		date > $com.hout
		cat $com | andromeda >> $com.hout
		date >> $com.hout
		sleep 2
		if ( $com =~ "*A.riz") then
			parse.csh $com
		endif
		if ( $com =~ "*Afound") then
			parse.csh found/A.riz
		endif
	end                                                                               
endif

