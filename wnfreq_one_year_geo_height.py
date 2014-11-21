"""
File: wnfreq_one_year_geo_height.py
Author: Charlie White
Email: charlie.white@mail.utoronto.ca
Github: echarliewhite
Description: This script is a test of the wnfreq_routines_2.0.py functions on
one year of ERA-Interim geopotential height data on one pressure surface.
Should be run on animus or sila.
"""
import wnfreq_routines_2_0 as wnfreq
import numpy as np
from matplotlib import pyplot as plt
import netCDF4 as nc4

g=9.81 # for conversion from geopotential to geopotential height

# import reanalysis data
filename = '/users/jk/14/cwhite/ERAInterim/dailymean/300hPa/z/z3d_2013_dailymean_NA_300hPa.nc'
var_name='Z_GDS0_ISBL'
nc = nc4.Dataset(filename)
data = nc.variables[var_name][:]
data /= g # convert from geopotential to geopotential height
lon = nc.variables['g0_lon_3'][:]
lat = nc.variables['g0_lat_2'][:]
plev = nc.variables['lv_ISBL1'][:]
time = nc.variables['initial_time0_hours'][:]
nc.close()

import pdb; pdb.set_trace()
