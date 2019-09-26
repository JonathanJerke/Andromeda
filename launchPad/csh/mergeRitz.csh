#!/bin/csh

#source launch.csh

if $#argv == 3 then
set src = $1
set dir = $2
set prev = $src/D.kry
set file = $dir/A.riz
echo "*Body $file" > $file
echo "*InputOutput" >> $file
echo "control found" >> $file
echo "read $dir/stage " >> $file
foreach f (`cat $3`)

#if ( `splitFlag.csh` ) then
echo "read $prev.$f-vector" >> $file
#else 
#echo "vector $prev.$f" >> $file
#endif

end
echo "control ritzPhase" >> $file
echo ".InputOutput" >> $file
echo "*Parameters" >> $file
echo "  collect 0" >> $file
echo ".Parameters" >> $file

echo ".Body" >> $file
${LAUNCH}/csh/go.csh $file &

endif
