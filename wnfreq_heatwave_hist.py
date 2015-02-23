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
g = 9.81 # for conversion from geopotential to geopotential height

def fft_1d_freq_wn(fft_coeffs_in,wn_max,N):
    """This function takes a 2D fft, computes the inverse fft, then computes 1D
    ffts along each axis, averaged over the non-fft axis. I.e. 2D fft in
    wavenumber and frequency becomes 2 1D ffts in wavenumber and frequency
    respectively.
    """
    real_space_geo = wnfreq.invert_wnfreq_spectrum(fft_coeffs_in,1,wn_max,N,tol=1e6)
    fft_coeffs_wn = np.fft.fft(real_space_geo,axis=1)
    fft_coeffs_wn = np.average(np.absolute(fft_coeffs_wn),axis=0)[1:wn_max+1]
    fft_coeffs_freq = np.fft.fft(real_space_geo,axis=0)
    fft_coeffs_freq = np.fft.fftshift(np.average(np.absolute(fft_coeffs_freq),axis=1))

    return fft_coeffs_freq, fft_coeffs_wn

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

# initialize arrays to hold FFTs for each N_day day period preceeding heat waves
# 3 arrays to hold average coefficients over wavenumbers
# 3 arrays to hold average coefficients over frequencies

hw_fft_total_wn = np.zeros(wn_max)
hw_fft_trav_wn = np.zeros(wn_max)
hw_fft_stand_wn = np.zeros(wn_max)
hw_fft_total_freq = np.zeros(N_days)
hw_fft_trav_freq = np.zeros(N_days)
hw_fft_stand_freq = np.zeros(N_days)

all_days_fft_total_wn = np.zeros(wn_max)
all_days_fft_trav_wn = np.zeros(wn_max)
all_days_fft_stand_wn = np.zeros(wn_max)
all_days_fft_total_freq = np.zeros(N_days)
all_days_fft_trav_freq = np.zeros(N_days)
all_days_fft_stand_freq = np.zeros(N_days)

Total_N_hw = 0 # track total number of heatwaves across all years
Total_N_ints = 0 # track total number of 20 day periods across all years

# for each year
for k in range(t_data.shape[0]):
    # determine which days have heat waves
    heatwave = fhw.find_heat_waves(t_data, 97.5, k, t_lat, t_lon, 44, 5)
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

    N_hw = np.sum(heatwavestart) # number of heatwaves
    N_ints = z_time.size/N_days # max possible number of FFT intervals

    # starting N_days into the time period, check for days where heatwaves start.
    # Take the N_days before that and take the FFT, store in arrays (3).
    for i in range(N_days,heatwavestart.size):
        if heatwavestart[i] == True:
            hw_fft_total_temp, hw_fft_stand_temp, hw_fft_trav_temp = \
              wnfreq.calc_wnfreq_spectrum(geoheight[i+shift_days-N_days:i+shift_days,:],wn_max)
            Total_N_hw += 1
            # calculate 1D freq/wn ffts
            fft_total_freq, fft_total_wn = fft_1d_freq_wn(hw_fft_total_temp,wn_max,z_lon.size)
            fft_stand_freq, fft_stand_wn = fft_1d_freq_wn(hw_fft_stand_temp,wn_max,z_lon.size)
            fft_trav_freq, fft_trav_wn = fft_1d_freq_wn(hw_fft_trav_temp,wn_max,z_lon.size)
            # add to total
            hw_fft_total_wn += fft_total_wn
            hw_fft_stand_wn += fft_stand_wn
            hw_fft_trav_wn += fft_trav_wn

            hw_fft_total_freq += fft_total_freq
            hw_fft_stand_freq += fft_stand_freq
            hw_fft_trav_freq += fft_trav_freq

    # now take same FFTs for all 20 day periods in the interval. (general climate)
    # consider taking all 20 day intervals NOT covered by heatwaves.
    for n in range(N_ints):
        i = n*N_days
        all_days_fft_total_temp, all_days_fft_stand_temp, all_days_fft_trav_temp = \
          wnfreq.calc_wnfreq_spectrum(geoheight[i:i+N_days,:],wn_max)
        Total_N_ints += 1
        # calculate 1D freq/wn ffts
        fft_total_freq, fft_total_wn = fft_1d_freq_wn(all_days_fft_total_temp,wn_max,z_lon.size)
        fft_stand_freq, fft_stand_wn = fft_1d_freq_wn(all_days_fft_stand_temp,wn_max,z_lon.size)
        fft_trav_freq, fft_trav_wn = fft_1d_freq_wn(all_days_fft_trav_temp,wn_max,z_lon.size)
        # add to total
        all_days_fft_total_wn += fft_total_wn
        all_days_fft_stand_wn += fft_stand_wn
        all_days_fft_trav_wn += fft_trav_wn

        all_days_fft_total_freq += fft_total_freq
        all_days_fft_stand_freq += fft_stand_freq
        all_days_fft_trav_freq += fft_trav_freq



    print float(k)*100/t_data.shape[0], '%'

# add up wavelength amplitudes for standing travelling total waves in
# pre-heatwave periods and divide by the number of periods (average amplitude
# at each wavelength). Take absolute value.
hw_fft_total_wn_average = hw_fft_total_wn/Total_N_hw
hw_fft_stand_wn_average = hw_fft_stand_wn/Total_N_hw
hw_fft_trav_wn_average  = hw_fft_trav_wn/Total_N_hw

hw_fft_total_freq_average = hw_fft_total_freq/Total_N_hw
hw_fft_stand_freq_average = hw_fft_stand_freq/Total_N_hw
hw_fft_trav_freq_average  = hw_fft_trav_freq/Total_N_hw

all_days_fft_total_wn_average = all_days_fft_total_wn/Total_N_ints
all_days_fft_stand_wn_average = all_days_fft_stand_wn/Total_N_ints
all_days_fft_trav_wn_average  = all_days_fft_trav_wn/Total_N_ints

all_days_fft_total_freq_average = all_days_fft_total_freq/Total_N_ints
all_days_fft_stand_freq_average = all_days_fft_stand_freq/Total_N_ints
all_days_fft_trav_freq_average  = all_days_fft_trav_freq/Total_N_ints

# look at start and -20 days? progression in 20 day period? average across 20 days?

# wavenumber FFT plot
bar_width = 0.35
index = np.arange(wn_max) + 1
fig,ax = plt.subplots(3,1)

ax[0].bar(index,all_days_fft_total_wn_average,width=bar_width,label='All Days')
ax[0].bar(index+bar_width,hw_fft_total_wn_average,width=bar_width,
          color='r',label='Pre-Heat Wave Days')
ax[0].set_title("Total Fourier Coeffs")
ax[0].legend()
ax[0].set_xticks(index,minor=True)
ax[1].bar(index,all_days_fft_stand_wn_average,width=bar_width,label='All Days')
ax[1].bar(index+bar_width,hw_fft_stand_wn_average,width=bar_width,
          color='r',label='Pre-Heat Wave Days')
ax[1].set_title("Standing Fourier Coeffs")
ax[1].legend()
ax[1].set_xticks(index,minor=True)
ax[2].bar(index,all_days_fft_trav_wn_average,width=bar_width,label='All Days')
ax[2].bar(index+bar_width,hw_fft_trav_wn_average,width=bar_width,
          color='r',label='Pre-Heat Wave Days')
ax[2].set_title("Travelling Fourier Coeffs")
ax[2].legend()
ax[2].set_xticks(index,minor=True)
ax[2].set_xlabel('Wavenumber')

# frequency FFT plot
fig,ax = plt.subplots(3,1)
index = np.arange(N_days) + 1

ax[0].bar(index,all_days_fft_total_freq_average,width=bar_width,label='All Days')
ax[0].bar(index+bar_width,hw_fft_total_freq_average,width=bar_width,
          color='r',label='Pre-Heat Wave Days')
ax[0].set_title("Total Fourier Coeffs")
ax[0].legend()
ax[0].set_xticks(index,minor=True)
ax[1].bar(index,all_days_fft_stand_freq_average,width=bar_width,label='All Days')
ax[1].bar(index+bar_width,hw_fft_stand_freq_average,width=bar_width,
          color='r',label='Pre-Heat Wave Days')
ax[1].set_title("Standing Fourier Coeffs")
ax[1].legend()
ax[1].set_xticks(index,minor=True)
ax[2].bar(index,all_days_fft_trav_freq_average,width=bar_width,label='All Days')
ax[2].bar(index+bar_width,hw_fft_trav_freq_average,width=bar_width,
          color='r',label='Pre-Heat Wave Days')
ax[2].set_title("Travelling Fourier Coeffs")
ax[2].legend()
ax[2].set_xticks(index,minor=True)
ax[2].set_xlabel('Frequency')
import pdb; pdb.set_trace()
plt.show()
# Other things to consider: change the period length, shift relative to start
# day, non-pre-heatwave days only, plot the geopotential at heatwave start
# dates/animate, different latitudes, take a histogram counting which
# wavenumber had max amplitude, eastward/westward
