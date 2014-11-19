#! /bin/bash

# This code to be used on animus/sila

dir="/climdata/ERAInterim/all_levels/dailymean/t/"
files=$dir"*nc"
# if directory changes change filename index
namestart=44

suffix="_summersurface.nc"

for file in $files
do
    filename=${file:$namestart}
    shortname=${filename%.nc}
    outputfile="/users/jk/14/cwhite/ERAInterim/dailymean/$shortname$suffix" 
    ncks -d t,[june],[august] -v t $file $outputfile
done
