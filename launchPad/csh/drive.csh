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

set flagBegin = 0
set flagEnd = 0
set begin = null
set end = null

if $#argv == 1 then
set flagBegin = 1
set begin = $1
endif

if $#argv == 2 then
set flagBegin = 1
set flagEnd = 0
set begin = $1
set end = $2
endif

if $#argv == 3 then
set flagBegin = 1
set flagEnd = 0
set begin = $1
set end = $2
set exclude = $3
endif


foreach com (`cat commands`)
        if ( $flagBegin && $com =~ "$begin*" ) then
                set flagBegin = 0
        endif

        if ( (! $flagEnd) && $com =~ "$end*" ) then
                set flagEnd = 1
        endif


        if ( (! $flagBegin )&& (! $flagEnd )&&(! ($com =~ "$exclude*"))) then              
	
                $LAUNCH/csh/go.csh $com
                echo "$com" >> driver
                sleep 3
                echo "$com"

        endif

end
