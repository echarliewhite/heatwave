#! /bin/bash

# This code to be used on animus/sila

files="/climdata/ERAInterim/all_levels/dailymean/t/*nc"

for file in $files
do
    filename=`expr match "$file" '.*\(t3d*nc\)'`
    outputfile="/users/jk/14/cwhite/ERAInterim/dailymean/$filename" 
    echo $outputfile
done
