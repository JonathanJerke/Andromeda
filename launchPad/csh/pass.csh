#!/bin/csh

set RS = RS
if ( -e $LAUNCH/references ) then
rm $LAUNCH/references
endif
foreach n  ( 1 2 )

foreach sa ( `cat $RS/sa-$n`)


foreach rs ( `cat $RS/rs`)
	if ($argv[1] == found ) then
		mkdir $RS-$n-$sa-$rs
		foreach p (`cat $RS/interactions/reg`)
                        mkdir $RS-$n-$sa-$rs/$p
                        ln  $RS/interactions/$p/* $RS-$n-$sa-$rs/$p
                end

		cp $RS/cnfg-$n $RS-$n-$sa-$rs/cnfg
 		cp $RS/boot-$n $RS-$n-$sa-$rs/boot
		cp $RS/cat-$n-$sa $RS-$n-$sa-$rs/cat
                cp $RS/post $RS-$n-$sa-$rs/post

	endif
	cp $RS/inc-$n $RS-$n-$sa-$rs/inc
 	foreach inc ( $argv )
		if ( $inc =~ :* ) then
 			set ii = $inc:s/://
 			cp $RS/$ii-$n $RS-$n-$sa-$rs/$ii
 		endif
	end
	cd $RS-$n-$sa-$rs
	if ($argv[1] == found ) then
			echo "*Body" 					>  wignerSeitz
			echo "*Parameters" 				>> wignerSeitz
			echo "	electronGasDensity $rs" 		>> wignerSeitz
			echo ".Parameters" 				>> wignerSeitz
			echo ".Body" 					>> wignerSeitz
			start.csh go $sa  
	endif
        mv *-* go
	
		cd go
			pwd >> $LAUNCH/references
			pathFinder.csh $argv   
		cd ..
	cd ..
end
end 
end
