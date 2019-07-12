#!/bin/csh

if ( $#argv == 2 )then

mkdir $1
chdir $1
mkdir found
cp ../../control/*found found
cp ../../control/stage found
cp ../../control/boot .
cp ../../symmetry/$2 symmetry
cp ../../control/states .
#$LAUNCH/csh/go.csh found/found found/Afound
else

	echo "run symmetry"
endif
