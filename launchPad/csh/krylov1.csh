#!/bin/csh

if $#argv == 4 then
set src = $1
set dir = $2
set prev = $src/D.kry
foreach i (`cat $3`)
foreach h (`cat $4`)
set file = $dir/B.kry.$h.$i
echo "*Body $file" > $file
echo "*InputOutput" >> $file
echo "read found/found" >> $file
echo "read $dir/stage " >> $file
echo "solo $prev.$i" >> $file
echo "control kry1Phase ">> $file
echo ".InputOutput" >> $file
echo "*Parameters" >> $file
echo "spam $h" >> $file
echo ".Parameters" >> $file
echo ".Body" >> $file
$LAUNCH/csh/go.csh $file &

end
sleep 1
end
else 
echo "originDir targetDir goSet"
endif

