"""
File: wnfreq_heatwave_hist.py
Author: Charlie White
Email: charlie.white@mail.utoronto.ca
Github: echarliewhite
Description: This script takes FFTs of geopotential height on identified
heatwave days and on regular days, then compares the distribution of
wavelengths using histograms.
"""
import wnfreq_routines_2_0 as wnfreq
import numpy as np
from matplotlib import pyplot as plt
import netCDF4 as nc4

g=9.81 # for conversion from geopotential to geopotential height

# import reanalysis data
filename = '/Users/charliewhite/Documents/Year4/Thesis/climdata/ERAInterim/dailymean/z/300hPa/z3d_2013_dailymean_NA_300hPa.nc'
var_name='Z_GDS0_ISBL'
nc = nc4.Dataset(filename)
data = nc.variables[var_name][:]
data /= g # convert from geopotential to geopotential height
lon = nc.variables['g0_lon_3'][:]
lat = nc.variables['g0_lat_2'][:]
plev = nc.variables['lv_ISBL1'][:]
time = nc.variables['initial_time0_hours'][:]
nc.close()

geoheight = np.squeeze(data) # remove lat and plev dimensions (single-valued)

wn_max = 20
