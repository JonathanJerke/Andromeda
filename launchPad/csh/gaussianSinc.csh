#!/bin/csh

if $#argv == 1 then
set src = $1
set dir = $1
set file = $dir/B.gaS
echo "*Body $file" > $file
echo "*InputOutput" >> $file
echo "read ../../control/found" >> $file
echo "read character ">> $file
echo "read $dir/stage " >> $file
echo "read ../../control/ritzPhase" >> $file
echo ".InputOutput" >> $file
echo ".Body" >> $file
${LAUNCH}/csh/go.csh $file &

endif
