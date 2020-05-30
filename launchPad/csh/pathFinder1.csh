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
#NOT USED				$LAUNCH/csh/hamiltonian.csh $curr
#NOT USED				$LAUNCH/csh/distill.csh $curr
#NOT USED				$LAUNCH/csh/distillExternal.csh $curr	
			sleep 3
				$LAUNCH/csh/gos.csh found/found
				$LAUNCH/csh/go.csh found/Afound
#				$LAUNCH/csh/prevKrylov.csh $curr $curr states
#		   		 sleep 3
#				$LAUNCH/csh/prevRitzB.csh $curr $curr states
	    		$LAUNCH/csh/prevKrylovC.csh  $curr $curr states
			endif 
		else
			$LAUNCH/csh/addStage.csh $prev $curr
			#define current iinc
			if ( -e $iinc ) then	
				cp $iinc $curr/inc
			else
				if ( -e ../$iinc ) then
					cp ../$iinc $curr/inc
				else
					exit
				endif
			endif
			#echo "cp $curr $iinc"
#NOT USED            		$LAUNCH/csh/hamiltonian.csh $curr
#NOT USED			$LAUNCH/csh/distill.csh $curr
#NOT USED                        $LAUNCH/csh/distillExternal.csh $curr
            		$LAUNCH/csh/krylov1.csh $prev $curr states spam
			sleep 3
			$LAUNCH/csh/sum1.csh $curr $curr states spam
            $LAUNCH/csh/ritz1.csh $prev $curr states
         	$LAUNCH/csh/prevKrylovB.csh  $curr $curr states
		endif
		set prev = $curr
		wait
	endif
end
wait
rm waiter
