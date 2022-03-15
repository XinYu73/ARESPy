#!/bin/bash 
# testing variables 
days=10 
guest="Katie"
echo "$guest checked in $days days ago" 
days=5 
guest="Jessica" 
echo "$guest checked in $days"

testing=$(date)
echo $testing
testing=$(who)
echo $testing

today=$(date +%m%d)
ls -al > log.$today

echo $[1+23]

var1=$(echo "scale=4; 3.44 / 5" | bc)
echo $var1

echo $(echo "$var1 * $var1" | bc)

var5=$(bc << EOF
scale=4
a1=($var1*$var1)
a1*a1
EOF
)

echo $var5
