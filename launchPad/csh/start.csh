#!/bin/csh

if ( $#argv == 2 ) then


mkdir $1
chdir $1
mkdir found
cp $LAUNCH/control/*found found
cp $LAUNCH/control/stage found

if ( -e ../boot ) then
cp ../boot .
else
cp $LAUNCH/control/boot .
endif

if ( -e ../inc ) then
cp ../inc .
else
cp $LAUNCH/control/inc .
endif

cp ../../symmetry/$2 symmetry

if ( -e ../post ) then
	cp ../post .
else 
	echo "none" > post
endif
echo "PostProcessing `cat post`"

set found = 10

if ( -e ../cnfg ) then
set body = `getParam.csh body ../cnfg`
set states = `getParam.csh states boot`
else
	exit
endif




if ( $body == 1 ) then
        switch ( $2 )
                case "all"
 			if ( -e cat ) then
			else
				genCat.csh 1 $found > ../cat
			endif
                        seq 1 1 > catalog
                        seq 1 $states  > states
                        breaksw
        endsw


endif


if ( $body == 2 ) then
        switch ( $2 )
                case "A1" 
                        if ( -e cat ) then
                        else
                                genCat.csh 4 $found > ../cat
                        endif

			seq 1 4 > catalog
                        seq 1 $states  > states
                        breaksw
                case "A2"
                        if ( -e cat ) then
                        else
                                genCat.csh 4 $found > ../cat
                        endif

                        seq 1 4> catalog
                        seq 1 $states  > states
                        breaksw
        endsw


endif



if ( $body == 3 ) then
	switch ( $2 ) 
		case "A1"

                        if ( -e cat ) then
                        else
                                genCat.csh 11 $found > ../cat
                        endif


			seq 1 11 > catalog
			seq 1 $states  > states
			breaksw
		case "A2"
                        if ( -e cat ) then
                        else
                                genCat.csh 11 $found > ../cat
                        endif

			seq 1 11 > catalog
			seq 1 $states > states
			breaksw
		case "E"
                        if ( -e cat ) then
                        else
                                genCat.csh 19 $found > ../cat
                        endif

			seq 1 19 > catalog
			seq 1 $states > states
			breaksw
	endsw


endif

if ( $body == 4 ) then
        switch ( $2 )
                case "A1"
                        if ( -e cat ) then
                        else
                                genCat.csh 39 $found > ../cat
                        endif

                        seq 1 39 > catalog
                        seq 1 $states  > states
                        breaksw
                case "A2"
                        if ( -e cat ) then
                        else
                                genCat.csh 39 $found > ../cat
                        endif

                        seq 1 39 > catalog
                        seq 1 $states > states
                        breaksw
                case "E"
                        if ( -e cat ) then
                        else
                                genCat.csh 55 $found > ../cat
                        endif

                        seq 1 55 > catalog
                        seq 1 $states  > states
                        breaksw
                case "T1"
                        if ( -e cat ) then
                        else
                                genCat.csh 58 $found > ../cat
                        endif

                        seq 1 58 > catalog
                        seq 1 $states  > states
                        breaksw
                case "T2"
                        if ( -e cat ) then
                        else
                                genCat.csh 58 $found > ../cat
                        endif

                        seq 1 58 > catalog
                        seq 1 $states  > states
                        breaksw


        endsw


endif


else

        echo "run symmetry"
endif
