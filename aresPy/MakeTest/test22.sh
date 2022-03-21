#!/bin/bash

#hompound comparisions
if [ -d $HOME ] && [ -w $HOME/testing ]
then
	echo "both condition satisfied"
elif [ -d $HOME ] || [ -w $HOME/testing ]
then
	echo "At least one is satisfied"
fi
if [ -d $HOME ] && [ -w $HOME/testing ]
then
	echo "both condition satisfied"
elif [ -d $HOME ] || [ -w $HOME/testing ]
then
	echo "At least one is satisfied"
fi

#
var1=10
if (( $var1 ** 2 > 90 ))
then
	(( var2 = $var1 ** 2 ))
	echo "The square of $var1 is $var2"
fi

if [[ $USER == l* ]]
then
	echo "Hello $USER"
fi

case $USER in
xinyu | luoxs)
	echo "well well"  
	echo "hello $USER";; # end a case block with ;;
testing)
	echo "For testing";;
*)
	echo "illegal";;
esac # end of case 
