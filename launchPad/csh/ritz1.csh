#!/bin/csh

#source launch.csh

if $#argv == 3 then
set src = $1
set dir = $2
set prev = $src/D.kry
set curr = $dir/b.sum
set file = $dir/C.riz
echo "*Body $file" > $file
echo "*InputOutput" >> $file
echo "read found/found" >> $file
echo "read $dir/stage " >> $file
foreach f (`cat $3`)
echo "vector $prev.$f" >> $file
echo "vector $curr.$f" >> $file
end
echo "control ritzPhase" >> $file
echo ".InputOutput" >> $file
echo "*Parameters" >>$file
echo "	collect 1 " >> $file
echo ".Parameters" >> $file
echo ".Body" >> $file
${LAUNCH}/csh/go.csh $file &

endif
