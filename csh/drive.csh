#!/bin/csh

foreach com ( `cat commands` )
		$LAUNCH/csh/go.csh $com
	sleep 3
	echo "$com"
end
