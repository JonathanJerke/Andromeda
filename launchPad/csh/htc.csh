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
$LAUNCHER_DIR/paramrun
