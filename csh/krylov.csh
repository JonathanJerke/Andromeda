#!/bin/csh

if $#argv == 3 then
set src = $1
set dir = $2
set prev = $src/A.riz
foreach i (`cat $3`)
set file = $dir/B.kry.${i}
echo "*Body $file" > $file
echo "*InputOutput" >> $file
echo "read ../../control/found" >> $file
echo "read $dir/stage " >> $file
#echo "vector $prev-${i}" >> $file
echo "read ../../control/kryPhase ">> $file
echo ".InputOutput" >> $file
echo ".Body" >> $file
$LAUNCH/csh/cats.csh $file $prev-$i &
end
else 
echo "originDir targetDir goSet"
endif

