#!/bin/sh


remote="$1"
url="$2"


while read local_ref local_sha remote_ref remote_sha
do
	if [[ $ref =~ .*/master$ ]];
	then
		echo "Run test before pushing to the master branch"
		exit 1
	fi
done

exit 0
