#!/bin/csh

if ( -e stop ) then
	echo "rm'd stop"
	rm stop
endif
if ( -e waiter ) then
	echo "rm'd waiter"
	rm waiter
endif
if ( -e driver ) then
	rm driver
endif

foreach com (`cat commands` )
	$LAUNCH/csh/go.csh $com
	echo "$com" >> driver
	sleep 3
	echo "$com"
end
