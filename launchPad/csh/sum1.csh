#!/bin/csh

if $#argv == 4 then
set src = $1
set dir = $2
set prev = $src/B.kry
foreach i (`cat $3`)
set file = $dir/b.sum.$i
echo "*Body $file" > $file
echo "*InputOutput" >> $file
echo "read found/found" >> $file
echo "read $dir/stage " >> $file
foreach h (`cat $4`)
echo "vector $prev.$h.$i" >> $file
end
echo "control sumPhase ">> $file
echo ".InputOutput" >> $file
echo ".Body" >> $file
$LAUNCH/csh/go.csh $file &

end
else 
echo "originDir targetDir goSet"
endif

