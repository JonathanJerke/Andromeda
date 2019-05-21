#!/bin/csh

#source launch.csh
#set LAUNCH = "/home/jjerke/launch"
if $#argv == 3 then
set dir = $2
cp character $dir ## keep a track in snow!
set src = $1
foreach i (`cat $3` )
echo $i
sleep 1
set prev = $src/A.riz
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

set prev2 = $dir/B.dcp
set file2 = $dir/C.kry.$i
echo "*Body $file2" > $file2
echo "*InputOutput" >> $file2
echo "read ../../control/found" >> $file2
echo "read character " >> $file2
echo "read $dir/stage " >> $file2
echo "vector $prev2.$i" >> $file2
echo "read ../../control/kryPhase ">> $file2
echo ".InputOutput" >> $file2
echo ".Body" >> $file2

${LAUNCH}/csh/go.csh $file $file2 &

#cat $file | olympics >> $file.hout  &
end
else
echo "originDir targetDir goSet"
endif
