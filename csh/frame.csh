#!/bin/csh

if $#argv == 1 then
set dir = $1
set src = $1
set prev = $src/D.fed.gs
set file = $dir/D.fed
rm *cell*
head -n 1 $1/B.kry.1.vector > $1/D.fed.gs.vector
echo "*Body $file" > $file
echo "*InputOutput" >> $file
echo "read ../../control/found" >> $file
#echo "read character ">> $file
echo "read $dir/stage " >> $file
echo "operator $prev" >> $file
endif
echo "read ../../control/fedPhase " >> $file
echo ".InputOutput" >> $file
echo ".Body" >> $file
${LAUNCH}/csh/go.csh $file &
else
echo "originDir targetDir goSet"
endif
