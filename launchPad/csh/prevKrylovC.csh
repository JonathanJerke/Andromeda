#!/bin/csh

if $#argv == 3 then
set src = $1
set dir = $2
set prev = $src/A.riz
foreach i (`cat $3`)
set file = $dir/D.kry.${i}
echo "*Body $file" > $file
echo "*InputOutput" >> $file
echo "read found/found" >> $file
echo "read $dir/stage " >> $file
#if ( `splitFlag.csh` ) then

#else
   echo "vector $prev-${i}" >> $file
#endif
echo "control buildPhase ">> $file
echo ".InputOutput" >> $file
echo ".Body" >> $file

#if ( `splitFlag.csh` ) then
#      		 $LAUNCH/csh/cats.csh $file $prev-$i 
#       else
		$LAUNCH/csh/go.csh $file	   
#  endif
end
else 
echo "originDir targetDir goSet"
endif

