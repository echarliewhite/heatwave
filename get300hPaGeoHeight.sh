#! /bin/bash

# This code to be used on animus/sila

dir="/climdata/ERAInterim/all_levels/dailymean/z/"
files=$dir"*nc"
# if directory changes change filename index
namestart=44

suffix="_NA_300hPa.nc"

for file in $files
do
    filename=${file:$namestart}
    shortname=${filename%.nc}
    outputfile="/users/jk/14/cwhite/ERAInterim/dailymean/300hPa/z/$shortname$suffix" 
    # may 1 is day 121, sept 30 is day 273
    # 300 hPa
    # whole circumference (along line of latitude)
    # latitude given by g0_lat_2
    ncks -d initial_time0_hours,121,273 -d g0_lat_2,40.0 -d lv_ISBL1,300.0 $file $outputfile
    echo $outputfile
done
