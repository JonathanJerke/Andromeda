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
	
	#		if ( `splitFlag.csh`) then
	#			$LAUNCH/csh/cats.csh found/found
	#		else
	#		endif
			$LAUNCH/csh/hamiltonian.csh $curr
			$LAUNCH/csh/distill.csh $curr
			sleep 3
			$LAUNCH/csh/gos.csh found/found
			$LAUNCH/csh/go.csh found/Afound
	#		if ( `splitFlag.csh`) then
     #                   $LAUNCH/csh/catKrylov.csh $curr $curr states
     #                   $LAUNCH/csh/mergeRitzB.csh $curr $curr states
	#		else
                       $LAUNCH/csh/prevKrylov.csh $curr $curr states
		       			sleep 3

                        $LAUNCH/csh/prevRitzB.csh $curr $curr states
	#		endif
	                $LAUNCH/csh/prevKrylovB.csh  $curr $curr states
      	
		endif 
	else
			$LAUNCH/csh/addStage.csh $prev $curr
                        $LAUNCH/csh/hamiltonian.csh $curr
			$LAUNCH/csh/distill.csh $curr
         #               if ( `splitFlag.csh`) then
          #              $LAUNCH/csh/splitKrylov.csh $prev $curr states
           #             $LAUNCH/csh/mergeRitzB.csh $curr $curr states
            #            else
                        $LAUNCH/csh/krylov.csh $prev $curr states
						sleep 3

                        $LAUNCH/csh/prevRitzB.csh $curr $curr states
             #           endif
                        $LAUNCH/csh/prevKrylovB.csh  $curr $curr states
	endif
	set prev = $curr
	#sleep 4
	wait
end
wait
#sort commands > commands.temp
#mv commands.temp commands
rm waiter
