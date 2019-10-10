#!/bin/csh

set flag = 1

foreach i (`grep catalog boot`)
	if ( $i =~ "*catalog*") then
	else	
		set flag = 0
		echo $i
	endif
end
if ( $flag ) then 
	echo "0"
endif
