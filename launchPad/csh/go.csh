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
		if ( `grep FINIS $com.hout` =~ "FINIS*" ) then
		else
			touch stop
			echo "stopped... now zombie"
		endif

		if ( $com =~ "*A.riz") then
	                if ( `splitFlag.csh`) then
				parse.csh $com
			else 
		#		copy.csh $com
			endif
		endif
		if ( $com =~ "*Afound") then
	              if ( `splitFlag.csh`) then
				parse.csh found/A.riz
		      else
		#	copy.csh found/A.riz
	              endif
		endif
		if ( $com =~ "*B.kry*") then
			sort -n -r $com.vector > $com.tm
			mv $com.tm $com.vector	
		endif
		sleep 2
	end                                                                               
endif

