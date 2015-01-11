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
import find_heat_waves as fhw
import numpy as np
from matplotlib import pyplot as plt
import netCDF4 as nc4

# User-defined variables
wn_max = 20 # max wave number to retain in FFT

# Constants
g=9.81 # for conversion from geopotential to geopotential height

# import reanalysis data
# pressure data
z_filename = '/Users/charliewhite/Documents/Year4/Thesis/climdata/ERAInterim/dailymean/z/300hPa/z3d_2013_dailymean_NA_300hPa.nc'
z_var_name='Z_GDS0_ISBL'
z_nc = nc4.Dataset(z_filename)
z_data = nc.variables[z_var_name][:]
z_data /= g # convert from geopotential to geopotential height
z_lon = z_nc.variables['g0_lon_3'][:]
z_lat = z_nc.variables['g0_lat_2'][:]
z_plev = z_nc.variables['lv_ISBL1'][:]
z_time = z_nc.variables['initial_time0_hours'][:]
z_nc.close()

geoheight = np.squeeze(data) # remove lat and plev dimensions (single-valued)

# temperature data
t_dir = '/Users/charliewhite/Documents/Year4/Thesis/climdata/ERAInterim/dailymean/t/summersurface/'
t_data, t_time, t_plev, t_lat, t_lon = fhw.import_temps(t_dir)

# determine which days have heat waves
heatwave = (t_data, 97.5, t_data.shape[0]-1) # 2013
# find the starting dates of heat wave events
heatwavestart = np.empty(heatwave.size,dtype=bool)
if heatwave[0] == True:
    heatwavestart[0] = True
else:
    heatwavestart = False
for i in range(1,heatwave.size):
    if heatwave[i] == True and heatwave[i-1] == False:
        heatwavestart[i] = True
    else:
        heatwavestart[i] = False

# initialize arrays (3) to hold FFTs for each 20 day period proceeding heat waves
# determine length of FFT array returned for each 30 day period
# for each 30 day period the array size is the same as the input array
# [time,...,lon], but the lon dimension is set to wn_max (higher wavenumbers
# are cut off).

# starting 20 days into the time period, check for days where heatwaves start.
# Take the 20 days before that and take the FFT, store in arrays (3).

# now take same FFTs for all 20 day periods in the interval. (general climate)
# consider taking all 20 day intervals NOT covered by heatwaves.

# add up wavelength amplitudes for standing travelling total waves in
# pre-heatwave periods and divide by the number of periods (average amplitude
# at each wavelength
# look at start and -20 days? progression in 20 day period? average across 20 days?

# do the same for the whole interval set.

# plot!

# Other things to consider: change the period length, shift relative to start day, non-pre-heatwave days only, plot the geopotential at heatwave start dates/animate, different latitudes
