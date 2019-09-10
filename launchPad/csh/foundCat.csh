#!/bin/csh

set printFlag = 0
foreach i (  $1:s/found/ /:s/found/ / ) #
	if ( $printFlag == 1 ) then
		echo $i:s/C/ /:s/C/ /
	else
		set printFlag = 1
	endif	
end

