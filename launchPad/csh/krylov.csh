#!/bin/csh

if $#argv == 3 then
set src = $1
set dir = $2
set prev = $src/D.kry
foreach i (`cat $3`)
echo "$i"
set file = $dir/B.kry.$i
echo "*Body $file" > $file
echo "*InputOutput" >> $file
echo "read found/found" >> $file
echo "read $dir/stage " >> $file
echo "vector $prev.$i" >> $file
echo "control kryPhase ">> $file
echo ".InputOutput" >> $file
echo ".Body" >> $file
$LAUNCH/csh/go.csh $file &

end
else 
echo "originDir targetDir goSet"
endif

