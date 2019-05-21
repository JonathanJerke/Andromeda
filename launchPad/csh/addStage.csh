#!/bin/csh

if $#argv == 5 then

set attack = $1
set stage = $2
set newRank = $3
set dir = $5
set newStage = $dir/stage
set src = $4
mkdir $dir
if ! -e $newStage  then

echo "*Body print " > $newStage
echo "*InputOutput" >> $newStage
echo " read $src/stage ">> $newStage
echo ".InputOutput" >> $newStage
echo "*Parameters " >> $newStage
echo "    attack        $attack " >> $newStage
echo "    bandStage     $stage " >> $newStage
echo "    basisStage    $newRank " >> $newStage
echo ".Parameters" >> $newStage
echo ".Body" >> $newStage

else
echo    "cannot overwrite targetDir"
endif

else
echo "attack bandStage basisStage originDir targetDir"
endif


