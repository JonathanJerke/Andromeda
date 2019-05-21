#!/bin/csh

#source launch.csh

if $#argv == 3 then
set src = $1
set dir = $2
set prev = $src/C.kry
set file = $dir/A.riz
echo "*Body $file" > $file
echo "*InputOutput" >> $file
echo "read ../../control/found" >> $file
echo "read character ">> $file
echo "read $dir/stage " >> $file
foreach f (`cat $3`)
echo "vector $prev.$f" >> $file
end
echo "read ../../control/ritzPhase" >> $file
echo ".InputOutput" >> $file
echo ".Body" >> $file
#cat $file | olympics
${LAUNCH}/csh/go.csh $file &

endif
