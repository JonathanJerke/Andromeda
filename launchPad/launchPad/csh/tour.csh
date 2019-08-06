#!/bin/csh

cd $LAUNCH

foreach dir (`cat references$1`)
	cd $dir
	drive.csh >> touring$1
	rm commands
	cd $LAUNCH
end
