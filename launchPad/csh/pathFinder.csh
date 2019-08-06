#!/bin/csh

if ( -e commands ) then
rm commands
endif

touch waiter
set flag = 1
foreach curr ( $argv )
	if ( $flag == 1 ) then
		set flag = 0
		if ( $curr == found ) then
	
			if ( `splitFlag.csh`) then
				$LAUNCH/csh/cats.csh found/found
			else
				$LAUNCH/csh/gos.csh found/found
			endif
			$LAUNCH/csh/go.csh found/Afound   
			$LAUNCH/csh/krylov.csh $curr $curr states
			$LAUNCH/csh/ritzB.csh $curr $curr states
	                $LAUNCH/csh/krylovB.csh  $curr $curr states
      	
		endif 
	else
		$LAUNCH/csh/addStage.csh $prev $curr
		$LAUNCH/csh/ritz.csh $prev $curr states
		$LAUNCH/csh/krylov.csh  $curr $curr states
		$LAUNCH/csh/ritzB.csh $curr $curr states
		$LAUNCH/csh/krylovB.csh  $curr $curr states
	endif
	set prev = $curr
	sleep 4
end
sleep 4
rm waiter
