#!/bin/csh

if ( $#argv == 1 ) then
set RS = $1
rm $LAUNCH/references2
foreach sp ( `cat $RS/sp`)
foreach cf ( `cat $RS/bt`)
foreach ic ( `cat $RS/ic`)
foreach sa ( `cat $RS/sa`)
set dir = $RS-$sa-$sp-$cf-$ic 
	mkdir $dir
	cp $RS/cnfg-$sp $dir/cnfg
 	cp $RS/boot-$cf $dir/boot
        cp $RS/inc-$ic $dir/inc
		cd $dir
		start.csh go $sa  
		cd go
			pwd >> $LAUNCH/references2
			pathFinder.csh found 1 2 3
		cd ..
	cd ..
end
end
end
end
else 
echo "RS dir"
endif
