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
		sleep 1
		if ( `grep FINIS $com.hout` =~ "FINIS*" ) then
		#	echo "finis"
			else
			echo "$com" > stop
			
		echo "stopped...on $com now zombie"
		endif
# allow D to be NORMAL!!
		if ( $com =~ "*.riz") then
	               # if ( `splitFlag.csh`) then
			#	echo "parse"
			#	parse.csh $com
			#else 
			#endif
		endif
		if ( $com =~ "*Afound") then
	          #    if ( `splitFlag.csh`) then
				parse.csh found/A.riz
		  #    else
		#	copy.csh found/A.riz
	          #    endif
		endif
		if ( $com =~ "*B.kry*" ) then
# && (`grep shift boot ` =~ "flow *" ||   `grep shift boot ` =~ "twist *"|| `grep shift boot ` =~ "shift *")) then
		sleep 1
			post.csh $com.vector > $com.tm
			mv $com.tm $com.vector
			sleep 1	
		endif
		sleep 2
	end                                                                               
endif

