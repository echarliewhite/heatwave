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
    outputfile="/users/jk/14/cwhite/ERAInterim/dailymean/t/summersurface/$shortname$suffix" 
    # may 1 is day 121, sept 30 is day 273
    # 1000 hPa
    # longitude from -125 e to -70 e, latitude from 25 n to 60 n (North America)
    ncks -d initial_time0_hours,121,273 -d g0_lon_3,-125.0,-70.0 -d g0_lat_2,25.0,50.0 -d lv_ISBL1,1000.0 $file $outputfile
    echo $outputfile
done
