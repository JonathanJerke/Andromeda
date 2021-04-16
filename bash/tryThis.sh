#!/bin/bash

#make sure the defs.sh are in order under a sort -n, I prefer A,B,C,D,E

pathFinder.sh found 1 2 3 4 5
cp commands  summa
summa.sh 5 6X
pathFinder.sh 6X 7 8 9 10 11 12 
grep E commands | tail -n 36 > summa
summa.sh 12 13X
pathFinder.sh 13X `seq 14 24`
grep E commands | tail -n 36 > summa
summa.sh 24 25X
pathFinder.sh 25X `seq 26 36`
cp inc1 7/inc
cp inc1 14/inc
cp inc1 21/inc
sort -n commands > c 
mv c commands
