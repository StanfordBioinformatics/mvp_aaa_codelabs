#!/usr/bin/bash

file=$1

for x in `cat $file` ; 
do s=`echo ${x##*/}` ; 
	echo $s ; { 
	time gsutil -m cp -L $s.gsutil.log $x/*/vcfs/* gs://gbsc-gcp-project-mvp-va_aaa/data/ ; 
} &> $s.log ; done
