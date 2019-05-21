#!/bin/csh

#source launch.csh
if $#argv == 1 then
set src = found
set dir = $1
set prev = $src/A.riz
set file = $dir/A.riz
echo "*Body $file" > $file
echo "*InputOutput" >> $file
echo "read ../../control/found" >> $file
echo "read character ">> $file
echo "read $dir/stage " >> $file
#foreach f (`cat goSet`)
#echo "vector $prev.$f" >> $file
#en
echo "vector found/found" >> $file
echo "read ../../control/ritzPhase" >> $file
echo ".InputOutput" >> $file
echo ".Body" >> $file
#cat $file | olympics
${LAUNCH}/csh/go.csh $file  &
else
echo "targetDir"

endif
