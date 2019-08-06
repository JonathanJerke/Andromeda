#!/bin/csh


set fileV = $1-vector

echo "*Body" > $fileV
echo "*InputOutput" >> $fileV
echo "	vector $1" >> $fileV
echo ".InputOutput" >> $fileV
echo ".Body" >> $fileV

$LAUNCH/csh/go.csh $1
