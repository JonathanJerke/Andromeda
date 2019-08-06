#!/bin/csh

if $#argv == 3 then
set src = $1
set dir = $2
set prev = $src/A.riz
foreach i (`cat $3`)
set file = $dir/C.prt.${i}
echo "*Body $file" > $file
echo "*InputOutput" >> $file
echo "read ../../control/found" >> $file
echo "read $dir/stage " >> $file
echo "vector $prev-${i}" >> $file
echo "read ../../control/kryPhase ">> $file
echo ".InputOutput" >> $file
echo "*Parameters" >> $file
echo "	filter 0 ">> $file
echo "	iterations 1 ">>$file
echo ".Parameters">>$file
echo ".Body" >> $file
$LAUNCH/csh/go.csh $file &
end
else
echo "originDir targetDir goSet"
endif

