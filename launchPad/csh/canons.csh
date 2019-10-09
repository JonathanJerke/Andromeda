#!/bin/csh

if $#argv == 2 then
set src = $1
set file = $src.ana.$2
echo "*Body $file" > $file
echo "*InputOutput" >> $file
echo "  read $src" >> $file
echo ".InputOutput" >> $file
echo "*Parameters ">> $file
echo "  iterations 1">> $file
echo "  decompose 4 ">> $file
echo "  basisRank  $2">> $file
echo ".Parameters" >> $file
echo ".Body" >> $file
$LAUNCH/csh/go.csh $file
else
echo "vector-file ##canoncial-ranks"
endif

