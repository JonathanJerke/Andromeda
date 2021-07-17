#!/bin/bash

for  com in $* 
	do

		set d1 = `head -n 1 $com`
		set t1 = `date -d "$d1" +%s `

		set d2 = `tail -n 1 $com`
		set t2 = `date -d "$d2" +%s `
		echo "$com $t1 $t2" | awk '{print $2 - $1}'

	done
