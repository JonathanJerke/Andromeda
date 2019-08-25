#!/bin/csh

if ( $#argv == 2 ) then

if ( -e boot && -e inc ) then 

if ( -e cnfg ) then
set body = `getParam.csh body cnfg`
set states = `getParam.csh states boot`
else
exit
endif

mkdir $1
chdir $1
mkdir found
cp ../../control/*found found
cp ../../control/stage found
cp ../boot .
cp ../inc .
cp ../../symmetry/$2 symmetry
echo "none" > post
echo "PostProcessing `cat post`"

if ( $body == 1 ) then
        switch ( $2 )
                case "all"
                        seq 1 1 > catalog
                        seq 1 $states  > states
                        breaksw
        endsw


endif


if ( $body == 2 ) then
        switch ( $2 )
                case "A1"
                        seq 1 4 > catalog
                        seq 1 $states  > states
                        breaksw
                case "A2"
                        seq 1 4> catalog
                        seq 1 $states  > states
                        breaksw
        endsw


endif



if ( $body == 3 ) then
	switch ( $2 ) 
		case "A1"
			seq 1 11 > catalog
			seq 1 $states  > states
			breaksw
		case "A2"
			seq 1 11 > catalog
			seq 1 $states > states
			breaksw
		case "E"
			seq 1 19 > catalog
			seq 1 $states > states
			breaksw
	endsw


endif

if ( $body == 4 ) then
        switch ( $2 )
                case "A1"
                        seq 1 39 > catalog
                        seq 1 $states  > states
                        breaksw
                case "A2"
                        seq 1 39 > catalog
                        seq 1 $states > states
                        breaksw
                case "E"
                        seq 1 55 > catalog
                        seq 1 $states  > states
                        breaksw
                case "T1"
                        seq 1 58 > catalog
                        seq 1 $states  > states
                        breaksw
                case "T2"
                        seq 1 58 > catalog
                        seq 1 $states  > states
                        breaksw


        endsw


endif
else
	
echo "place boot and inc with cnfg"
endif



else

        echo "run symmetry"
endif
