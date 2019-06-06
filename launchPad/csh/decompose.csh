#!/bin/csh

#source launch.csh
#set LAUNCH = "/home/jjerke/launch"
if $#argv == 3 then
set dir = $2
set src = $1
set prev = $src/A.riz
foreach i (`cat $3` )
echo $i
sleep 1
set file = $dir/B.dcp.$i
echo "*Body $file" > $file
echo "*InputOutput" >> $file
echo "read ../../control/found" >> $file
echo "read character ">> $file
echo "read $dir/stage " >> $file
echo "vector $prev-$i" >> $file
echo "read ../../control/dcpPhase " >> $file
echo ".InputOutput" >> $file
echo ".Body" >> $file

${LAUNCH}/csh/go.csh $file &

#cat $file | olympics >> $file.hout  &
end
else
echo "originDir targetDir goSet"
endif