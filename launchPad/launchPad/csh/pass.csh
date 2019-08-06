#!/bin/csh


rm $LAUNCH/references
foreach n  ( 2 )

foreach sa ( `cat sa-$n`)


foreach rs ( `cat rs`)
	mkdir gas-$n-$sa-$rs
	cp cnfg-$n gas-$n-$sa-$rs/cnfg
		cd gas-$n-$sa-$rs
		echo "*Body" 					>  wignerSeitz
		echo "*Parameters" 				>> wignerSeitz
		echo "	electronGasDensity $rs" >> wignerSeitz
		echo ".Parameters" 				>> wignerSeitz
		echo ".Body" 					>> wignerSeitz
		start.csh go $sa $n 2
		cd go
			pwd >> $LAUNCH/references
			pathFinder.csh found 1 2 3 4 
		cd ..
	cd ..
end
end 
end
