#!/bin/csh
module load launcher
setenv LAUNCHER_JOB_FILE launcher
setenv LAUNCHER_NHOSTS 1
setenv LAUNCHER_WORKDIR `pwd`
if (-e ppn) then
	setenv LAUNCHER_PPN `cat ppn`
else
	echo "set PPN for speedup"
	setenv LAUNCHER_PPN 4
endif

if ( -e jobDone ) then
rm jobDone
endif

echo "touch jobDone" >> launcher

foreach i in ( `seq 1 100`)

        $LAUNCHER_DIR/paramrun
        if ( -e jobDone ) then
                exit
        else
                sleep 100
        endif
end
