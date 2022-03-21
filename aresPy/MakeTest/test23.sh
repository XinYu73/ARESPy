#!/bin/bash
for name in $(ls)
do
	echo "There is $name"
done

list="Alabama Alaska Arizona Arkansas Colorado"
list=$list" Connectichut"  #concatenate
for state in $list
do
	echo "We got $state"
done

IFSOLD=$IFS
IFS=$'\n'
for word in $(cat test22.sh)
do
	echo "$word "
done
IFS=$IFSOLD

#reading a directory using wildcards
for file in ./*
do
	if [ -d $file ]
	then
		echo "$file is a directory"
	elif [ -f $file ]
	then
		echo "$file is a file"
	fi
done

