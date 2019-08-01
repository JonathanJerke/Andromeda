#!/bin/csh

if ( $#argv == 4 ) then

mkdir $1
chdir $1
mkdir found
cp ../../control/*found found
cp ../../control/stage found
cp ../../control/boot .
cp ../../symmetry/$2 symmetry


if ( $3 == 1 ) then
        switch ( $2 )
                case "all"
                        seq 1 1 > catalog
                        seq 1 $4  > states
                        breaksw
        endsw


endif


if ( $3 == 2 ) then
        switch ( $2 )
                case "A1"
                        seq 1 4 > catalog
                        seq 1 $4  > states
                        breaksw
                case "A2"
                        seq 1 4> catalog
                        seq 1 $4  > states
                        breaksw
        endsw


endif



if ( $3 == 3 ) then
	switch ( $2 ) 
		case "A1"
			seq 1 11 > catalog
			seq 1 $4  > states
			breaksw
		case "A2"
			seq 1 11 > catalog
			seq 1 $4 > states
			breaksw
		case "E"
			seq 1 19 > catalog
			seq 1 $4 > states
			breaksw
	endsw


endif

if ( $3 == 4 ) then
        switch ( $2 )
                case "A1"
                        seq 1 39 > catalog
                        seq 1 $4  > states
                        breaksw
                case "A2"
                        seq 1 39 > catalog
                        seq 1 $4 > states
                        breaksw
                case "E"
                        seq 1 55 > catalog
                        seq 1 $4  > states
                        breaksw
                case "T1"
                        seq 1 58 > catalog
                        seq 1 $4  > states
                        breaksw
                case "T2"
                        seq 1 58 > catalog
                        seq 1 $4  > states
                        breaksw


        endsw


endif

else

        echo "run symmetry BODY #states"
endif
