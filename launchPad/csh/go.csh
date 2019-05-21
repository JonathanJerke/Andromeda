#!/bin/csh

#source launch.csh
#if $#argv == 1 then

foreach com ( $argv )
date > $com.hout
cat $com | ~/Perci/olympics >> $com.hout 
sleep 2
date >> $com.hout
end
