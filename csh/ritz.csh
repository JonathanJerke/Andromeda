#!/bin/csh

#source launch.csh

if $#argv == 3 then
set src = $1
set dir = $2
set prev = $src/B.kry
set file = $dir/A.riz
echo "*Body $file" > $file
echo "*InputOutput" >> $file
echo "read ../../control/found" >> $file
#echo "read go ">> $file
echo "read $dir/stage " >> $file
foreach f (`cat $3`)
echo "read $prev.$f-vector" >> $file
end
echo "read ../../control/ritzPhase" >> $file
echo ".InputOutput" >> $file
echo ".Body" >> $file
#cat $file | olympics
${LAUNCH}/csh/go.csh $file &

endif
