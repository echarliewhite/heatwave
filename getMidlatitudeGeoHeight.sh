#! /bin/bash

# This code to be used on animus/sila

dir="/climdata/ERAInterim/all_levels/dailymean/z/"
files=$dir"*nc"
# if directory changes change filename index
namestart=44

suffix="_midlatitude.nc"

for file in $files
do
    filename=${file:$namestart}
    shortname=${filename%.nc}
    outputfile="/users/jk/14/cwhite/ERAInterim/dailymean/z/NH/$shortname$suffix" 
    # may 1 is day 121, sept 30 is day 273
    # whole circumference (along lines of latitude)
    # latitude given by g0_lat_2
    # whole troposphere + lower stratosphere included
    ncks -d initial_time0_hours,121,273 -d g0_lat_2,0.0,90.0 -d lv_ISBL1,50.0,1000.0 $file $outputfile
    echo $outputfile
done
