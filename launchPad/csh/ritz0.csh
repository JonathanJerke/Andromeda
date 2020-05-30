#!/bin/csh

#source launch.csh

if $#argv == 4 then
set src = $1
set dir = $2
set prev = $src/D.kry
set curr = $dir/B.kry
set file = $dir/C.riz
echo "*Body $file" > $file
echo "*InputOutput" >> $file
echo "read found/found" >> $file
echo "read $dir/stage " >> $file
foreach f (`cat $3`)
foreach h (`cat $4`)
echo "vector $curr.$h.$f" >> $file
end
echo "vector $prev.$f" >> $file
end
echo "control ritzPhase" >> $file
echo ".InputOutput" >> $file
echo "*Parameters" >>$file
echo "	collect 1 " >> $file
echo ".Parameters" >> $file
echo ".Body" >> $file
${LAUNCH}/csh/go.csh $file &

endif
