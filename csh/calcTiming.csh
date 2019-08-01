#!/bin/csh

if $#argv == 1 then

set d1 = `head -n 1 $1`
set t1 = `date -d "$d1" +%s `

set d2 = `tail -n 1 $1`
set t2 = `date -d "$d2" +%s `
echo "$t1 $t2" | awk '{print $2 - $1}'

endif
