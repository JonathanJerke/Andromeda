#!/bin/csh

#set LAUNCH = `which LAUNCH`
#alaunch.csh
if $#argv == 3 then
set src = $1
set dir = $2
set prev = $src/B.dcp
foreach i (`cat $3`)
echo "$i"
sleep 1
set file = $dir/C.kry.$i
echo "*Body $file" > $file
echo "*InputOutput" >> $file
echo "read ../../control/found" >> $file
echo "read character " >> $file
echo "read $dir/stage " >> $file
echo "vector $prev.$i" >> $file
echo "read ../../control/kryPhase ">> $file
echo ".InputOutput" >> $file
echo ".Body" >> $file
#cat $file | olympics & 
$LAUNCH/csh/go.csh $file &

end
else 
echo "originDir targetDir goSet"
endif
