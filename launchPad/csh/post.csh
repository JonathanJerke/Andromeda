#!/bin/csh

if ( -e post ) then
	if ( `cat post` =~ "none" ) then
		cat $1 
	endif
	if ( `cat post` =~ "reverse" ) then
		sort -n -r $1 
	endif
	if ( `cat post` =~ "last" ) then
		tail -n 1 $1
	endif
	if ( `cat post` =~ "lastFew" ) then
                tail -n 3 $1
        endif
	
	if ( `cat post` =~ "copyLast" ) then

                tail -n 1 $1 > $1:s/B/D/     
		tail -n 1 $1                                                            
	endif      	

else
	echo "oops set post"
	echo "post" > stop
endif

