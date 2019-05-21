#!/bin/csh

if $#argv == 6 then

set shift = $1
set attack = $2
set stage = $3
set newRank = $4
set dir = $6
set newStage = $dir/stage
set src = $5
mkdir $dir
if ! -e $newStage  then

echo "*Body print " > $newStage
echo "*InputOutput" >> $newStage
echo " read $src/stage ">> $newStage
echo ".InputOutput" >> $newStage
echo "*Parameters " >> $newStage
echo "    shift 	$shift	" >> $newStage
echo "    attack        $attack " >> $newStage
echo "    bandStage     $stage " >> $newStage
echo "    basisStage    $newRank " >> $newStage
echo ".Parameters" >> $newStage
echo ".Body" >> $newStage

else
echo    "cannot overwrite targetDir"
endif

else
echo "shift attack bandStage basisStage originDir targetDir"
endif


