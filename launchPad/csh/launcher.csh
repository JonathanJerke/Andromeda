#!/bin/csh

if ( -e stop ) then
        echo "rm'd stop"
        rm stop
endif
if ( -e waiter ) then
        echo "rm'd waiter"
        rm waiter
endif
if ( -e launcher ) then
        rm launcher
endif

set com0 = `head -n 1 commands`
set h = `header.csh $com0`
foreach com (`cat commands`)
	set ch = `header.csh $com`
	if ( $h =~ $ch ) then

	else 
		echo `wc -l launcher ` 
		if ( `wc -l launcher ` =~ "1 *" ) then
       			set COM = `cat launcher`
			$COM
			echo "$COM"
		else 
			echo "htc $h"
			htc.csh
		endif
		rm launcher
		set h = $ch
	endif	             
	echo "go.csh $com " >> launcher
end


                if ( `wc -l launcher ` == 1 ) then
                        foreach COM ( `cat launcher` )
                                go.csh $COM
                                echo "go $COM"
                        end
                else 
                        echo "htc $h"
                        htc.csh
                endif
                rm launcher

