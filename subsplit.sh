#!/bin/bash

f=$1
n=$2

while read -r line; do
	working=$(echo $line | grep -v ">")
	for y in `seq 0 ${#working}`; do
		stringZ=${working:y:n}
       		if [ ${#stringZ} == $n ]
		then 
			echo $stringZ
		fi
	done
done < $f
