#!/bin/csh

set printFlag = 0
foreach i (  $1:s/./ / ) 
	if ( $printFlag == 0 ) then
		echo $i
		set printFlag = 1
	endif	
end

