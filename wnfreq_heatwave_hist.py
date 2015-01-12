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
N_days = 20 # number of days over which to take FFT
shift_days = 0 # number of days to shift the FFT interval forward (0 -> ends
               # day before heatwave)

# Constants
g=9.81 # for conversion from geopotential to geopotential height

# import reanalysis data
# pressure data
z_filename = '/Users/charliewhite/Documents/Year4/Thesis/climdata/ERAInterim/dailymean/z/300hPa/z3d_2013_dailymean_NA_300hPa.nc'
z_var_name='Z_GDS0_ISBL'
z_nc = nc4.Dataset(z_filename)
z_data = z_nc.variables[z_var_name][:]
z_data /= g # convert from geopotential to geopotential height
z_lon = z_nc.variables['g0_lon_3'][:]
z_lat = z_nc.variables['g0_lat_2'][:]
z_plev = z_nc.variables['lv_ISBL1'][:]
z_time = z_nc.variables['initial_time0_hours'][:]
z_nc.close()

geoheight = np.squeeze(z_data) # remove lat and plev dimensions (single-valued)

# temperature data
t_dir = '/Users/charliewhite/Documents/Year4/Thesis/climdata/ERAInterim/dailymean/t/summersurface/'
t_data, t_time, t_plev, t_lat, t_lon = fhw.import_temps(t_dir)

# determine which days have heat waves
heatwave = fhw.find_heat_waves(t_data, 97.5, t_data.shape[0]-1) # 2013
# find the starting dates of heat wave events
heatwavestart = np.empty(heatwave.size,dtype=bool)
if heatwave[0] == True:
    heatwavestart[0] = True
else:
    heatwavestart[0] = False
for i in range(1,heatwave.size):
    if heatwave[i] == True and heatwave[i-1] == False:
        heatwavestart[i] = True
    else:
        heatwavestart[i] = False

# initialize arrays (3) to hold FFTs for each N_day day period proceeding heat waves
# for each N_day period the array size is the same as the input array
# [time,...,lon], but the lon dimension is set to wn_max (higher wavenumbers
# are cut off). Dimensions are therefore [N_hw,Ndays,wn_max]
N_hw = np.sum(heatwavestart) # number of heatwaves
N_ints = z_time.size/N_days # max possible number of FFT intervals
hw_fft_total = np.zeros((N_hw,N_days,wn_max))
hw_fft_trav = np.zeros((N_hw,N_days,wn_max))
hw_fft_stand = np.zeros((N_hw,N_days,wn_max))

all_days_fft_total = np.zeros((N_ints,N_days,wn_max))
all_days_fft_trav = np.zeros((N_ints,N_days,wn_max))
all_days_fft_stand = np.zeros((N_ints,N_days,wn_max))

# starting N_days into the time period, check for days where heatwaves start.
# Take the N_days before that and take the FFT, store in arrays (3).
n = 0 # for keeping track of the heatwave count
for i in range(20,heatwavestart.size):
    if heatwavestart[i] == True:
        hw_fft_total[n,:,:], hw_fft_stand[n,:,:], hw_fft_trav[n,:,:] = \
          wnfreq.calc_wnfreq_spectrum(geoheight[i+shift_days-N_days:i+shift_days,:],wn_max)
        n += 1

# now take same FFTs for all 20 day periods in the interval. (general climate)
# consider taking all 20 day intervals NOT covered by heatwaves.
for n in range(N_ints):
    i = n*N_days
    all_days_fft_total[n,:,:], all_days_fft_stand[n,:,:], all_days_fft_trav[n,:,:] = \
      wnfreq.calc_wnfreq_spectrum(geoheight[i:i+N_days,:],wn_max)

# add up wavelength amplitudes for standing travelling total waves in
# pre-heatwave periods and divide by the number of periods (average amplitude
# at each wavelength)
hw_fft_total_average = np.average(np.average(hw_fft_total,axis=0),axis=0)
hw_fft_stand_average = np.average(np.average(hw_fft_stand,axis=0),axis=0)
hw_fft_trav_average = np.average(np.average(hw_fft_trav,axis=0),axis=0)

all_days_fft_total_average = np.average(np.average(all_days_fft_total,axis=0),axis=0)
all_days_fft_stand_average = np.average(np.average(all_days_fft_stand,axis=0),axis=0)
all_days_fft_trav_average = np.average(np.average(all_days_fft_trav,axis=0),axis=0)

diff_fft_total_average = hw_fft_total_average - all_days_fft_total_average
diff_fft_stand_average = hw_fft_stand_average - all_days_fft_stand_average
diff_fft_trav_average = hw_fft_trav_average - all_days_fft_trav_average

# look at start and -20 days? progression in 20 day period? average across 20 days?

# plot
fig,ax = plt.subplots(3,3)
ax[0,0].plot(hw_fft_total_average)
ax[0,0].set_title("Pre-Heatwave Total FFT")
ax[0,1].plot(hw_fft_stand_average)
ax[0,1].set_title("Pre-Heatwave Standing FFT")
ax[0,2].plot(hw_fft_trav_average)
ax[0,2].set_title("Pre-Heatwave Travelling FFT")
ax[1,0].plot(all_days_fft_total_average)
ax[1,0].set_title("All Days Total FFT")
ax[1,1].plot(all_days_fft_stand_average)
ax[1,1].set_title("All Days Standing FFT")
ax[1,2].plot(all_days_fft_trav_average)
ax[1,2].set_title("All Days Travelling FFT")
ax[2,0].plot(diff_fft_total_average)
ax[2,0].set_title("Difference Total FFT")
ax[2,1].plot(diff_fft_stand_average)
ax[2,1].set_title("Difference Standing FFT")
ax[2,2].plot(diff_fft_trav_average)
ax[2,2].set_title("Difference Travelling FFT")

import pdb; pdb.set_trace()

plt.show()
# Other things to consider: change the period length, shift relative to start
# day, non-pre-heatwave days only, plot the geopotential at heatwave start
# dates/animate, different latitudes, take a histogram counting which
# wavenumber had max amplitude
