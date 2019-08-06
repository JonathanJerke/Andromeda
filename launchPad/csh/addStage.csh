#!/bin/csh

if $#argv == 2 then
set dir = $2
set newStage = $dir/stage
set src = $1
mkdir $dir
if ! -e $newStage  then

echo "*Body print " > $newStage
echo "*InputOutput" >> $newStage
echo " read $src/stage ">> $newStage
echo " read $dir/inc" >> $newStage
cp inc $dir
echo ".InputOutput" >> $newStage
echo ".Body" >> $newStage

else
echo    "cannot overwrite targetDir"
endif

else
echo "attack originDir targetDir"
endif


