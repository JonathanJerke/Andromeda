#!/bin/csh

if ( -e commands ) then
rm commands
endif

touch waiter
set flag = 1
foreach i ($argv)
        if ( ! $flag ) then
        $LAUNCH/csh/canons.csh $1 $i
        endif
        set flag = 0
end
wait
rm waiter

