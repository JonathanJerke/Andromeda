#!/bin/csh

if ( $#argv == 1 ) then
set RS = $1
rm $LAUNCH/references
foreach n  ( 3 4 )

foreach sa ( `cat $RS/sa-$n`)


foreach rs ( `cat $RS/rs`)
	mkdir $RS-$n-$sa-$rs
	cp $RS/cnfg-$n $RS-$n-$sa-$rs/cnfg
 	cp $RS/boot-$n $RS-$n-$sa-$rs/boot
        cp $RS/inc-$n $RS-$n-$sa-$rs/inc
	cp $RS/cat-$n-$sa $RS-$n-$sa-$rs/cat
		cd $RS-$n-$sa-$rs
		echo "*Body" 					>  wignerSeitz
		echo "*Parameters" 				>> wignerSeitz
		echo "	electronGasDensity $rs" 		>> wignerSeitz
		echo ".Parameters" 				>> wignerSeitz
		echo ".Body" 					>> wignerSeitz
		start.csh go $sa  
		cd go
			pwd >> $LAUNCH/references
			pathFinder.csh found 1   
		cd ..
	cd ..
end
end 
end
else 
echo "RS dir"
endif
