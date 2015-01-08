"""
File: wnfreq_one_year_geo_height.py
Author: Charlie White
Email: charlie.white@mail.utoronto.ca
Github: echarliewhite
Description: This script is a test of the wnfreq_routines_2.0.py functions on
one year of ERA-Interim geopotential height data on one pressure surface.
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

# average over latitudes and remove pressure dimension
datalon = np.squeeze(np.average(data,axis=2))

wn_max = 5 # max wavenumber for spectrum calculation
plot_freq_cutoff = 2
vertical_scale = 1

datatrans, datastand, dataprop = wnfreq.calc_wnfreq_spectrum(datalon,wn_max)

# plot FFT coeffs
wnfreq.plot_wnfreq_spectrum_lineplots(datatrans,1,plot_freq_cutoff,
                                      vertical_scale,leg_label='Total')
# plot standing FFT coeffs
wnfreq.plot_wnfreq_spectrum_lineplots(datastand,2,plot_freq_cutoff,
                               vertical_scale,leg_label='Standing')
# plot travelling FFT coeffs
wnfreq.plot_wnfreq_spectrum_lineplots(dataprop,3,plot_freq_cutoff,
                               vertical_scale,leg_label='Travelling')
plt.show()
