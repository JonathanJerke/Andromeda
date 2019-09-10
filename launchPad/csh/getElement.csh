#!/bin/csh

@ ct = 0
foreach i ( `cat $1`) 
	@ ct += 1
	#echo $ct
	if ( $ct == $2 ) then
		echo $i
	endif
end
