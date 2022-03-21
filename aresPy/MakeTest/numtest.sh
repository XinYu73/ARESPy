#!/bin/bash
value1=10
value2=12
if [ $value2 -ge $value1 ]
then
	echo "value2 is greater"
fi

testDir=/home/luoxs
if [ -d $testDir ]
then 
	echo "The $testDir directory exists"
else
	echo "$testDir doesn't exist"
fi

