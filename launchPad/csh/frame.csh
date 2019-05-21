#!/bin/csh

if $#argv == 3 then
set dir = $2
set src = $1
set prev = $src/B.dcp
set file = $dir/D.fed
echo "*Body $file" > $file
echo "*InputOutput" >> $file
echo "read ../../control/found" >> $file
echo "read character ">> $file
echo "read $dir/stage " >> $file

if 0 then
#have your own vector files to define density
else 

foreach i (`cat $3` )
echo "operator $prev.$i" >> $file
end

endif
echo "read ../../control/fedPhase " >> $file
echo ".InputOutput" >> $file
echo ".Body" >> $file
${LAUNCH}/csh/go.csh $file &
else
echo "originDir targetDir goSet"
endif
