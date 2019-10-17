#!/bin/csh

if ( -e commands ) then
rm commands
endif

touch waiter
set flag = 1
set iinc = inc
foreach curr ( $argv )
	if ( $curr =~ :* ) then
		set iinc = $curr:s/://
		#echo $iinc
	else 
		if ( $flag == 1 ) then
			set flag = 0
			if ( $curr == found ) then
				$LAUNCH/csh/hamiltonian.csh $curr
				$LAUNCH/csh/distill.csh $curr
				sleep 3
				$LAUNCH/csh/gos.csh found/found
				$LAUNCH/csh/go.csh found/Afound
				$LAUNCH/csh/prevKrylov.csh $curr $curr states
		   		 sleep 3
				$LAUNCH/csh/prevRitzB.csh $curr $curr states
	    		$LAUNCH/csh/prevKrylovB.csh  $curr $curr states
			endif 
		else
			$LAUNCH/csh/addStage.csh $prev $curr
			#define current iinc
			cp ../$iinc $curr/inc
			#echo "cp $curr $iinc"
            $LAUNCH/csh/hamiltonian.csh $curr
			$LAUNCH/csh/distill.csh $curr
            $LAUNCH/csh/krylov.csh $prev $curr states
			sleep 3
            $LAUNCH/csh/prevRitzB.csh $curr $curr states
         	$LAUNCH/csh/prevKrylovB.csh  $curr $curr states
		endif
		set prev = $curr
		wait
	endif
end
wait
rm waiter
