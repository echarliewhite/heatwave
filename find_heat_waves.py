"""
File: find_heat_waves.py
Author: Charlie White
Email: charlie.white@mail.utoronto.ca
Github: echarliewhite
Description: This script will identify heat waves based on a threshold
temperature defined relative to historical temperatures.
"""
import numpy as np
import netCDF4 as nc4
import os

g=9.81 # for conversion from geopotential to geopotential height

for file in os.listdir('/users/jk/14/cwhite/ERAInterim/dailymean/summersurface/t'):
    var_name='T_GDS0_ISBL'
    nc = nc4.Dataset(file)
    data = nc.variables[var_name][:]
    data /= g # convert from geopotential to geopotential height
    lon = nc.variables['g0_lon_3'][:]
    lat = nc.variables['g0_lat_2'][:]
    plev = nc.variables['lv_ISBL1'][:]
    time = nc.variables['initial_time0_hours'][:]
    nc.close()

import pdb; pdb.set_trace()
