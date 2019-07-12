#!/bin/csh

set att = 0.5
set bandStage = 1
set basisStage = 3

if ( -e commands ) then
rm commands
endif

touch waiter
set flag = 1
foreach curr ( $argv )
	if ( $flag == 1 ) then
		set flag = 0
		if ( $curr == found ) then
			$LAUNCH/csh/go.csh found/found found/Afound   
			$LAUNCH/csh/krylov.csh found found states          	
		endif 
	else
		$LAUNCH/csh/addStage.csh $att $bandStage $basisStage $prev $curr
		$LAUNCH/csh/ritz.csh $prev $curr states
		$LAUNCH/csh/krylov.csh  $curr $curr states		
	endif
	set prev = $curr
end
sleep 2
rm waiter
