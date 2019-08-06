#!/bin/csh

rm -rf found 
mkdir found
cp ../../control/*found found
$LAUNCH/csh/go.csh found/found found/Afound

