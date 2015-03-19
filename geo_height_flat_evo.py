"""
File: geo_height_flat_evo.py
Author: Charlie White
Email: charlie.white@mail.utoronto.ca
Github: echarliewhite
Description: This script produces charts similar to fig. 1 in
Kushner/Watt-Meyer decomposition paper for heat waves.
"""
import numpy as np
import heat_wave_tools as hwt
import wnfreq_routines_2_0 as wnfreq
import matplotlib.pyplot as plt
import sys
import gc

# import geoheight and temp data
z_data, z_time, z_plev, z_lat, z_lon = hwt.import_nc_dir(
    '/Users/charliewhite/Documents/Year4/Thesis/climdata/ERAInterim/dailymean/z/NH/',
    'Z_GDS0_ISBL', pbound=(300.0,None), latbound=(51.,24.))

g = 9.80665 # global average gravity
z_data /= g # convert from geopotential into geopotential height

t_data, t_time, t_plev, t_lat, t_lon = hwt.import_nc_dir(
    '/Users/charliewhite/Documents/Year4/Thesis/climdata/ERAInterim/dailymean/t/summersurface/',
    'T_GDS0_ISBL', pbound=(1000.0,None))

# get heat wave data
N_points = 0.05*t_lat.size*t_lon.size # 5% of grid points
heat_wave_dict = hwt.find_heat_waves(t_data,97.5,t_lat,t_lon,N_points,5)
print '%d heat waves detected' % len(heat_wave_dict.keys())

z_clim = np.average(z_data, axis=0) # climatology
# geoheight anomaly
z_anom = z_data
for i in range(z_data.shape[0]):
    z_anom[i] -= z_clim

# remove zonal mean (wavenumber-0 discarded by FFT)
z_anom_zonal = np.average(z_anom, axis=4)
for i in range(z_data.shape[4]):
    z_anom[:,:,:,:,i] -= z_anom_zonal

del z_data
gc.collect()

# Compute FFT and Inverse FFT
wn_max = 15
sys.stdout.write('Computing FFT...   ')
sys.stdout.flush()
z_anom_trans, z_anom_trans_standing, z_anom_trans_travelling = \
        wnfreq.calc_wnfreq_spectrum(z_anom, wn_max)
sys.stdout.write('Done\n')
del z_anom_trans
gc.collect()

sys.stdout.write('Computing Inverse FFT of Standing Component...   ')
sys.stdout.flush()
z_anom_standing = wnfreq.invert_wnfreq_spectrum(
        z_anom_trans_standing,1,wn_max,z_lon.size,tol=1e6)
sys.stdout.write('Done\n')
del z_anom_trans_standing
gc.collect()

sys.stdout.write('Computing Inverse FFT of Travelling Component...   ')
sys.stdout.flush()
z_anom_travelling = wnfreq.invert_wnfreq_spectrum(
        z_anom_trans_travelling,1,wn_max,z_lon.size,tol=1e6)
sys.stdout.write('Done\n')
del z_anom_trans_travelling
gc.collect()

print 'Preparing for plotting...'
# create composite geopotential for heat waves
comp_length = 60 # number of days centered on day 0 (even # only)

z_anom_comp = np.zeros((1,comp_length,1,z_lat.size,z_lon.size))
z_anom_comp_standing = np.zeros((1,comp_length,1,z_lat.size,z_lon.size))
z_anom_comp_travelling = np.zeros((1,comp_length,1,z_lat.size,z_lon.size))
center = np.zeros(2)
N_hw = 0


# find average heat wave center point
for hwave in heat_wave_dict.values():
    day0 = hwave.start
    # only accept fully captured heat waves
    if comp_length/2 <= day0 <= z_time.size - comp_length/2:
        center += np.array([hwave.lat,hwave.lon])
        N_hw += 1

center /= N_hw
center_ind = np.zeros(2,dtype=int)
z_lat_sort = np.argsort(z_lat)
z_lon_sort = np.argsort(z_lon)
# this function returns indices as if array sorted in ascending order,
# z_lat is sorted in descending order, so flip index about midpoint
center_ind[0] = (np.searchsorted(z_lat,center[0],sorter=z_lat_sort) - z_lat.size)*-1
center_ind[1] = np.searchsorted(z_lon,center[1],sorter=z_lon_sort)

print 'Shifting and averaging...'
# for plotting limits
inner_lat_bound = [0,z_lat.size]
inner_lon_bound = [0,z_lon.size]
for hwave in heat_wave_dict.values():
    day0 = hwave.start
    # only accept fully captured heat waves
    if comp_length/2 <= day0 <= z_time.size - comp_length/2:
        year = hwave.year
        hwave_lat = hwave.lat
        hwave_lon = hwave.lon
        
        # shift data to align at average center point
        lat_ind = (np.searchsorted(z_lat,hwave_lat,sorter=z_lat_sort) - z_lat.size)*-1
        lon_ind = np.searchsorted(z_lon,hwave_lon,sorter=z_lon_sort)
        lat_pad = lat_ind - center_ind[0]
        lon_pad = lon_ind - center_ind[1]
        lat_bound = [max(0,lat_pad),min(z_lat.size,z_lat.size+lat_pad)]
        lon_bound = [max(0,lon_pad),min(z_lon.size,z_lon.size+lon_pad)]
        inner_lat_bound = [max(lat_bound[0],inner_lat_bound[0]),
                           min(lat_bound[1],inner_lat_bound[1])]
        inner_lon_bound = [max(lon_bound[0],inner_lon_bound[0]),
                           min(lon_bound[1],inner_lon_bound[1])]
        z_anom_shift = np.zeros(z_anom_comp.shape)
        z_anom_shift_standing = np.zeros(z_anom_comp.shape)
        z_anom_shift_travelling = np.zeros(z_anom_comp.shape)
        z_anom_shift[:,:,:,lat_bound[0]:lat_bound[1],lon_bound[0]:lon_bound[1]] = \
            z_anom[year,day0-comp_length/2:day0+comp_length/2,:,
                    lat_bound[0]:lat_bound[1],lon_bound[0]:lon_bound[1]]
        z_anom_shift_standing[:,:,:,lat_bound[0]:lat_bound[1],lon_bound[0]:lon_bound[1]] = \
            z_anom_standing[year,day0-comp_length/2:day0+comp_length/2,:,
                    lat_bound[0]:lat_bound[1],lon_bound[0]:lon_bound[1]]
        z_anom_shift_travelling[:,:,:,lat_bound[0]:lat_bound[1],lon_bound[0]:lon_bound[1]] = \
            z_anom_travelling[year,day0-comp_length/2:day0+comp_length/2,:,
                    lat_bound[0]:lat_bound[1],lon_bound[0]:lon_bound[1]]

        # take absolute value
        z_anom_comp += np.absolute(z_anom_shift)
        z_anom_comp_standing += np.absolute(z_anom_shift_standing)
        z_anom_comp_travelling += np.absolute(z_anom_shift_travelling)

del z_anom
del z_anom_standing
del z_anom_travelling
del z_anom_shift
del z_anom_shift_standing
del z_anom_shift_travelling
gc.collect()
        
print 'Plotting...'
z_anom_comp /= N_hw
z_anom_comp_standing /= N_hw
z_anom_comp_travelling /= N_hw

# get limits on region of interest for amplitude
# midlatitude North America
amp_lat_bound = (np.array([np.searchsorted(z_lat,70,sorter=z_lat_sort),
                 np.searchsorted(z_lat,10,sorter=z_lat_sort)]) - z_lat.size)*-1
amp_lon_bound = [np.searchsorted(z_lon,-180,sorter=z_lon_sort),
                 np.searchsorted(z_lon,-30,sorter=z_lon_sort)]
import pdb; pdb.set_trace()


# average and max at each day
daily_mean_amp = np.zeros(z_anom_comp.shape[1])
daily_mean_amp_standing = np.zeros(z_anom_comp.shape[1])
daily_mean_amp_travelling = np.zeros(z_anom_comp.shape[1])
daily_peak_amp = np.zeros(z_anom_comp.shape[1])
daily_peak_amp_standing = np.zeros(z_anom_comp.shape[1])
daily_peak_amp_travelling = np.zeros(z_anom_comp.shape[1])
for i in range(z_anom_comp.shape[1]):
    daily_mean_amp[i] = np.mean(z_anom_comp[0,i,:,
        amp_lat_bound[0]:amp_lat_bound[1],amp_lon_bound[0]:amp_lon_bound[1]])
    daily_mean_amp_standing[i] = np.mean(z_anom_comp_standing[0,i,:,
        amp_lat_bound[0]:amp_lat_bound[1],amp_lon_bound[0]:amp_lon_bound[1]])
    daily_mean_amp_travelling[i] = np.mean(z_anom_comp_travelling[0,i,:,
        amp_lat_bound[0]:amp_lat_bound[1],amp_lon_bound[0]:amp_lon_bound[1]])
    daily_peak_amp[i] = np.max(z_anom_comp[0,i,:,
        amp_lat_bound[0]:amp_lat_bound[1],amp_lon_bound[0]:amp_lon_bound[1]])
    daily_peak_amp_standing[i] = np.max(z_anom_comp_standing[0,i,:,
        amp_lat_bound[0]:amp_lat_bound[1],amp_lon_bound[0]:amp_lon_bound[1]])
    daily_peak_amp_travelling[i] = np.max(z_anom_comp_travelling[0,i,:,
        amp_lat_bound[0]:amp_lat_bound[1],amp_lon_bound[0]:amp_lon_bound[1]])

days = np.arange(z_anom_comp.shape[1]) - z_anom_comp.shape[1]/2
plt.plot(days,daily_mean_amp,'r',label='Total Average Amplitude')
plt.plot(days,daily_mean_amp_standing,'b',label='Standing Average Amplitude')
plt.plot(days,daily_mean_amp_travelling,'g',label='Travelling Average Amplitude')
plt.plot(days,daily_peak_amp,'r--',label='Total Peak Amplitude')
plt.plot(days,daily_peak_amp_standing,'b--',label='Standing Peak Amplitude')
plt.plot(days,daily_peak_amp_travelling,'g--',label='Travelling Peak Amplitude')
plt.xlabel('Days +/- Heat Wave Start')
plt.ylabel('Geoheight Anomaly Amplitude [m]')
plt.title('Geoheight Anomaly Amplitude Near Heat Waves')
plt.legend()

plt.show()

