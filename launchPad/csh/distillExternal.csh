#!/bin/csh

#source launch.csh

if $#argv == 1 then
set dir = $1
set file = $dir/A.dix
echo "*Body $file" > $file
echo "*InputOutput" >> $file
echo "control found" >> $file
echo "read $dir/stage " >> $file
echo "control distillExternalPhase" >> $file
echo ".InputOutput" >> $file

echo ".Body" >> $file
${LAUNCH}/csh/go.csh $file &

endif