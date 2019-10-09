#!/bin/csh

set i = $1
set rl = `grep "383 " $i`
set r = $rl[2]:s/(//:s/)//

set dl = `grep omplete $i`
set d = $dl[3]

set hl = `grep State $i`
set h = $hl[2]:s/,//:s/,//

echo "$r  $d    $h "

