#!/bin/csh


set fileV = $1-vector

echo "*Body" > $fileV
echo "*InputOutput" >> $fileV




foreach c (`cat catalog`)
	set file = $1C${c}C
	echo "*Body $1C${c}C" > $file
	echo "*InputOutput" >> $file
	echo "read $1" >> $file
        if ( $#argv == 2 ) then
                echo "  vector $2C${c}C " >> $file
     	 endif

	echo ".InputOutput" >> $file
	echo "*Parameters" >> $file
	echo "	catalog $c" >> $file
	echo ".Parameters">> $file
	echo ".Body" >> $file
	
	echo "	vector $1C${c}C" >> $fileV
	
	$LAUNCH/csh/go.csh $1C${c}C
end
	
echo ".InputOutput" >> $fileV
echo ".Body" >> $fileV




