#!/bin/csh

foreach com ( `cat commands` )
	$LAUNCH/csh/go.csh $com
end
