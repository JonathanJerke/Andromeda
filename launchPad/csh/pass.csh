#!/bin/csh

if ( $#argv == 1 ) then
set RS = $1
rm $LAUNCH/references
foreach n  ( 2 3 )

foreach sa ( `cat $RS/sa-$n`)


foreach rs ( `cat $RS/rs`)
	mkdir $RS-$n-$sa-$rs
	cp $RS/cnfg-$n $RS-$n-$sa-$rs/cnfg
 	cp $RS/boot-$n $RS-$n-$sa-$rs/boot
        cp $RS/inc-$n $RS-$n-$sa-$rs/inc

		cd $RS-$n-$sa-$rs
		echo "*Body" 					>  wignerSeitz
		echo "*Parameters" 				>> wignerSeitz
		echo "	electronGasDensity $rs" 		>> wignerSeitz
		echo ".Parameters" 				>> wignerSeitz
		echo ".Body" 					>> wignerSeitz
		start.csh go $sa $n 6
		cd go
			pwd >> $LAUNCH/references
			pathFinder.csh found 1 2 3 
		cd ..
	cd ..
end
end 
end
else 
echo "RS dir"
endif
