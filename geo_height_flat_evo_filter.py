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

# climatology
z_clim = np.average(z_data, axis=0) 
t_clim = np.average(t_data, axis=0)

# anomaly
z_anom = z_data
t_anom = t_data
for i in range(z_data.shape[0]):
    z_anom[i] -= z_clim
for i in range(t_data.shape[0]):
    t_anom[i] -= t_clim

# remove zonal mean
z_anom_zonal = np.average(z_anom, axis=4)
t_anom_zonal = np.average(t_anom, axis=4)
for i in range(z_data.shape[4]):
    z_anom[:,:,:,:,i] -= z_anom_zonal
for i in range(t_data.shape[4]):
    t_anom[:,:,:,:,i] -= t_anom_zonal

del z_data
del t_data
del z_anom_zonal
del t_anom_zonal
gc.collect()

# Compute FFT and Inverse FFT
wn_max = 15
filter_wn = 5
z_filter = np.zeros((z_anom.shape[1],z_anom.shape[2],z_anom.shape[3],wn_max))
z_filter[:,:,:,filter_wn-1] = np.ones(z_filter.shape[0:3])
z_anom_standing = np.zeros(z_anom.shape, dtype='complex128')
z_anom_travelling = np.zeros(z_anom.shape, dtype='complex128')

for i in range(z_anom.shape[0]):
    sys.stdout.write('Computing FFT '+str(i+1)+'/'+str(z_anom.shape[0])+'...   ')
    sys.stdout.flush()
    z_anom_trans, z_anom_trans_standing, z_anom_trans_travelling = \
            wnfreq.calc_wnfreq_spectrum(z_anom[i], wn_max)

    z_anom_trans *= z_filter
    z_anom[i] = wnfreq.invert_wnfreq_spectrum(
            z_anom_trans,1,wn_max,z_lon.size,tol=1e6)
    sys.stdout.write('Done\n')
    del z_anom_trans
    gc.collect()

    sys.stdout.write('Computing Inverse FFT of Standing Component '+str(i+1)+'/'+str(z_anom.shape[0])+'...   ')
    sys.stdout.flush()
    z_anom_trans_standing *= z_filter
    z_anom_standing[i] = wnfreq.invert_wnfreq_spectrum(
            z_anom_trans_standing,1,wn_max,z_lon.size,tol=1e6)
    sys.stdout.write('Done\n')
    del z_anom_trans_standing
    gc.collect()

    sys.stdout.write('Computing Inverse FFT of Travelling Component '+str(i+1)+'/'+str(z_anom.shape[0])+'...   ')
    sys.stdout.flush()
    z_anom_trans_travelling *= z_filter
    z_anom_travelling[i] = wnfreq.invert_wnfreq_spectrum(
            z_anom_trans_travelling,1,wn_max,z_lon.size,tol=1e6)
    sys.stdout.write('Done\n')
    del z_anom_trans_travelling
    gc.collect()

print 'Preparing for plotting...'

comp_length = 60 # number of days centered on day 0
z_lat_sort = np.argsort(z_lat)
z_lon_sort = np.argsort(z_lon)
plot_lon_bound = [np.searchsorted(z_lon,-180,sorter=z_lon_sort),
                 np.searchsorted(z_lon,-30,sorter=z_lon_sort)]

# plot each heat wave
show = True # for checking format on first iteration
for hwave in heat_wave_dict.values():
    day0 = hwave.start
    # only accept fully captured heat waves
    if comp_length//2 <= day0 <= z_time.size - comp_length//2:
        year = hwave.year
        hwave_lon = hwave.lon

        z_plot_data = np.zeros((3,)+z_anom.shape,dtype='complex128')
        z_plot_data[0] = z_anom
        z_plot_data[1] = z_anom_standing
        z_plot_data[2] = z_anom_travelling
        fig, ax = hwt.plot_meridional_mean_evo(
                    z_plot_data[:,:,:,:,:,plot_lon_bound[0]:plot_lon_bound[1]],
                    z_lon[plot_lon_bound[0]:plot_lon_bound[1]], t_anom, t_lon,
                    year, day0, center_lon=hwave_lon,
                    day_range=[-comp_length//2,comp_length//2])
        fig.set_size_inches(17,11,forward=True)
        fig.suptitle('Year '+str(year)+' Day '+str(day0)+' Geopotential Height Anomaly')
        #if show==True: fig.show()
        #show==False
        fig.savefig('output/wn'+str(filter_wn)+'year'+str(year)+'day'+str(day0)+'.png',
                    format='png')

# plot composite
z_anom_comp = np.zeros((1,comp_length,1,z_lat.size,z_lon.size),dtype='complex128')
z_anom_comp_standing = np.zeros((1,comp_length,1,z_lat.size,z_lon.size),dtype='complex128')
z_anom_comp_travelling = np.zeros((1,comp_length,1,z_lat.size,z_lon.size),dtype='complex128')

t_anom_comp = np.zeros((1,comp_length,1,t_lat.size,t_lon.size),dtype='complex128')
t_anom_comp_standing = np.zeros((1,comp_length,1,t_lat.size,t_lon.size),dtype='complex128')
t_anom_comp_travelling = np.zeros((1,comp_length,1,t_lat.size,t_lon.size),dtype='complex128')

center = np.zeros(2)
N_hw = 0


# find average heat wave center point
for hwave in heat_wave_dict.values():
    day0 = hwave.start
    # only accept fully captured heat waves
    if comp_length//2 <= day0 <= z_time.size - comp_length//2:
        center += np.array([hwave.lat,hwave.lon])
        N_hw += 1

center /= N_hw
center_ind = np.zeros(2,dtype=int)
# this function returns indices as if array sorted in ascending order,
# z_lat is sorted in descending order, so flip index about midpoint
center_ind[0] = (np.searchsorted(z_lat,center[0],sorter=z_lat_sort) - z_lat.size)*-1
center_ind[1] = np.searchsorted(z_lon,center[1],sorter=z_lon_sort)

print 'Shifting and averaging...'
# for plotting limits
for hwave in heat_wave_dict.values():
    day0 = hwave.start
    # only accept fully captured heat waves
    if comp_length//2 <= day0 <= z_time.size - comp_length//2:
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
        t_lat_bound = lat_bound - np.array([0,z_lat.size - t_lat.size])
        t_lon_bound = lon_bound - np.array([0,z_lon.size - t_lon.size])

        z_anom_shift = np.zeros(z_anom_comp.shape,dtype='complex128')
        z_anom_shift_standing = np.zeros(z_anom_comp.shape,dtype='complex128')
        z_anom_shift_travelling = np.zeros(z_anom_comp.shape,dtype='complex128')
        z_anom_shift[:,:,:,lat_bound[0]:lat_bound[1],lon_bound[0]:lon_bound[1]] = \
            z_anom[year,day0-comp_length//2:day0+comp_length//2,:,
                    lat_bound[0]:lat_bound[1],lon_bound[0]:lon_bound[1]]
        z_anom_shift_standing[:,:,:,lat_bound[0]:lat_bound[1],lon_bound[0]:lon_bound[1]] = \
            z_anom_standing[year,day0-comp_length//2:day0+comp_length//2,:,
                    lat_bound[0]:lat_bound[1],lon_bound[0]:lon_bound[1]]
        z_anom_shift_travelling[:,:,:,lat_bound[0]:lat_bound[1],lon_bound[0]:lon_bound[1]] = \
            z_anom_travelling[year,day0-comp_length//2:day0+comp_length//2,:,
                    lat_bound[0]:lat_bound[1],lon_bound[0]:lon_bound[1]]

        t_anom_shift = np.zeros(t_anom_comp.shape)
        t_anom_shift[:,:,:,t_lat_bound[0]:t_lat_bound[1],t_lon_bound[0]:t_lon_bound[1]] = \
            t_anom[year,day0-comp_length//2:day0+comp_length//2,:,
                    t_lat_bound[0]:t_lat_bound[1],t_lon_bound[0]:t_lon_bound[1]]

        z_anom_comp += z_anom_shift
        z_anom_comp_standing += z_anom_shift_standing
        z_anom_comp_travelling += z_anom_shift_travelling

        t_anom_comp += t_anom_shift

z_anom_comp /= N_hw
z_anom_comp_standing /= N_hw
z_anom_comp_travelling /= N_hw

t_anom_comp /= N_hw

del z_anom
del z_anom_standing
del z_anom_travelling
del z_anom_shift
del z_anom_shift_standing
del z_anom_shift_travelling
gc.collect()
        
# plot
print 'Plotting...'
plot_lon_bound = [np.searchsorted(z_lon,-180,sorter=z_lon_sort),
                 np.searchsorted(z_lon,-30,sorter=z_lon_sort)]
z_plot_data = np.zeros((3,)+z_anom_comp.shape,dtype='complex128')
z_plot_data[0] = z_anom_comp
z_plot_data[1] = z_anom_comp_standing
z_plot_data[2] = z_anom_comp_travelling
fig, ax = hwt.plot_meridional_mean_evo(
                z_plot_data[:,:,:,:,:,plot_lon_bound[0]:plot_lon_bound[1]],
                z_lon[plot_lon_bound[0]:plot_lon_bound[1]], t_anom_comp, t_lon,
                0, comp_length//2, center_lon=center[1],
                day_range=[-comp_length//2,comp_length//2])
fig.savefig('output/march24composite.png',
            format='png')

plt.show()

