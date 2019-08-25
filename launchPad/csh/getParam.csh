#!/bin/csh

set string = `grep $1 $2`
shift string
echo $string
