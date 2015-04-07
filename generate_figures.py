"""
File: generate_figures.py
Author: Charlie White
Email: charlie.white@mail.utoronto.ca
Github: echarliewhite
Description: This script will generate all figures for my final document with
appropriate formatting/labels.
"""
import numpy as np
import heat_wave_tools as hwt
import wnfreq_routines_2_0 as wnfreq
import matplotlib.pyplot as plt
import sys
import gc
import datetime

pressure = 300.0 # the pressure level for all geoheight plots
# specific heat waves to plot together
plot_hwaves = ['year_9_day_49', 'year_32_day_32'] 
init_date = datetime.date(1979,5,1) # May 1, 1979

# import geoheight and temp data
z_data, z_time, z_plev, z_lat, z_lon = hwt.import_nc_dir(
    '/Users/charliewhite/Documents/Year4/Thesis/climdata/ERAInterim/dailymean/z/NH/',
    'Z_GDS0_ISBL', pbound=(pressure,None), latbound=(61.5,13.5))

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
del z_clim
del t_clim
del z_anom_zonal
del t_anom_zonal
gc.collect()

# plot boundary values
z_lat_sort = np.argsort(z_lat)
z_lon_sort = np.argsort(z_lon)
plot_lon_bound = [np.searchsorted(z_lon,-160.5,sorter=z_lon_sort),
                  np.searchsorted(z_lon,-28.5,sorter=z_lon_sort)]
"""Temperature Evolution Plots"""
# 3 days
plot_data = np.zeros((len(plot_hwaves),)+t_anom[:1,:5,:,:,:].shape)
date = np.empty(len(plot_hwaves), dtype=datetime.date)
center = np.zeros((len(plot_hwaves),2))
for i in range(len(plot_hwaves)):
    hwave = heat_wave_dict[plot_hwaves[i]]
    plot_data[i][0] = t_anom[hwave.year,hwave.start:hwave.start+5,:,:,:]
    center[i] = np.array([hwave.lat,hwave.lon])
    date[i] = datetime.date(init_date.year+hwave.year, init_date.month,
                            init_date.day) + datetime.timedelta(days=hwave.start)

fig, ax, cax, cbar = hwt.plot_evo(plot_data, 0, 0, t_lat, t_lon, center=center,
                                    days=[0,2,4], cmap='rainbow',
                                    clev_bound=[-15,15], clev_num=31,
                                    lat_ticks=5, lon_ticks=31)
plt.subplots_adjust(wspace=0.025)
datestring = date[0].isoformat()
ax[0,0].set_title('Heatwave starting on '+datestring, fontsize=22)
datestring = date[1].isoformat()
ax[0,1].set_title('Heatwave starting on '+datestring, fontsize=22)
fig.text(0.45, 0.06, 'Longitude', fontsize=18)
fig.text(0.05, 0.5, 'Latitude', rotation='vertical', fontsize=18)
cbar.set_label(r'Temperature Anomaly [$^\circ$C]', fontsize=18)

fig.savefig('figures/example_temperature_maps.png', format='png')
""""""

"""Geoheight Anomaly Maps"""
plot_data = np.zeros((len(plot_hwaves),) +\
              z_anom[:1, :13, :, :, plot_lon_bound[0]:plot_lon_bound[1]].shape)
for i in range(len(plot_hwaves)):
    hwave = heat_wave_dict[plot_hwaves[i]]
    plot_data[i][0] = z_anom[hwave.year, hwave.start-6:hwave.start+7, :, :,
                             plot_lon_bound[0]:plot_lon_bound[1]]

fig, ax, cax, cbar = hwt.plot_evo(plot_data[:1], 0, 6, z_lat,
                                z_lon[plot_lon_bound[0]:plot_lon_bound[1]],
                                center=center[:1], days=[0], show_day=False,
                                clev_bound=[-450,450], clev_num=37)
fig.set_size_inches(17,4.5)
datestring = date[0].isoformat()
ax.set_title(datestring, fontsize=22)
fig.text(0.45, 0.06, 'Longitude', fontsize=18)
fig.text(0.05, 0.5, 'Latitude', rotation='vertical', fontsize=18)
cbar.set_label(r'Geopotential Height Anomaly [m]', fontsize=14)

fig.savefig('figures/example_geoheight_map1.png', format='png')

fig, ax, cax, cbar = hwt.plot_evo(plot_data[1:], 0, 6, z_lat,
                                    z_lon[plot_lon_bound[0]:plot_lon_bound[1]],
                                    center=center[1:], days=[0], show_day=False,
                                    clev_bound=[-450,450], clev_num=37)
fig.set_size_inches(17,4.5)
datestring = date[1].isoformat()
ax.set_title(datestring, fontsize=22)
fig.text(0.45, 0.06, 'Longitude', fontsize=14)
fig.text(0.05, 0.5, 'Latitude', rotation='vertical', fontsize=18)
cbar.set_label(r'Geopotential Height Anomaly [m]', fontsize=14)

fig.savefig('figures/example_geoheight_map2.png', format='png')

fig, ax, cax, cbar = hwt.plot_evo(plot_data, 0, 6, z_lat,
                                    z_lon[plot_lon_bound[0]:plot_lon_bound[1]],
                                    center=center, days=[-6,-3,0,3,6],
                                    clev_bound=[-450,450], clev_num=37)
plt.subplots_adjust(wspace=0.025)
datestring = date[0].isoformat()
ax[0,0].set_title(datestring, fontsize=22)
datestring = date[1].isoformat()
ax[0,1].set_title(datestring, fontsize=22)
fig.text(0.45, 0.06, 'Longitude', fontsize=18)
fig.text(0.05, 0.5, 'Latitude', rotation='vertical', fontsize=18)
cbar.set_label(r'Geopotential Height Anomaly [m]', fontsize=18)

fig.savefig('figures/example_geoheight_maps.png', format='png')

""""""

comp_length = 60 # number of days centered on day 0
"""Geoheight Anomaly Meridional Average Plots"""
plot_data = np.zeros((len(plot_hwaves),) +\
              z_anom[:1, :comp_length, :, :, plot_lon_bound[0]:plot_lon_bound[1]].shape)

for i in range(len(plot_hwaves)):
    hwave = heat_wave_dict[plot_hwaves[i]]
    plot_data[i][0] = z_anom[hwave.year,
            hwave.start-comp_length//2:hwave.start+comp_length//2, :, :,
            plot_lon_bound[0]:plot_lon_bound[1]]
    
fig, ax, cax, cbar = hwt.plot_custom_meridional_mean_evo(plot_data,
        z_lon[plot_lon_bound[0]:plot_lon_bound[1]], 0, comp_length//2,
        center_lon=center[:,1], day_range=[-comp_length//2,comp_length//2],
        clev_bound1=[-260,260], clev_num1=27)

fig.set_size_inches(8.5,11)
datestring = date[0].isoformat()
ax[0].set_title(datestring, fontsize=22)
datestring = date[1].isoformat()
ax[1].set_title(datestring, fontsize=22)
fig.text(0.4, 0.055, 'Longitude', fontsize=18)
ax[0].set_ylabel('Days +/- Heat Wave Start', fontsize=18)
cbar.set_label(r'Meridional Mean Geopotential Height Anomaly [m]', fontsize=18)

fig.savefig('figures/example_geoheight_evolution_mean.png', format='png')
""""""

# Compute FFT and Inverse FFT
wn_max = 15
z_anom_standing = np.zeros(z_anom.shape, dtype='complex128')
z_anom_travelling = np.zeros(z_anom.shape, dtype='complex128')
# for filtering out 
filter_wn = 5
z_filter = np.zeros((z_anom.shape[1],z_anom.shape[2],z_anom.shape[3],wn_max))
z_filter[:,:,:,filter_wn-1] = np.ones(z_filter.shape[0:3])
z_anom_filter = np.zeros(z_anom.shape, dtype='complex128')
z_anom_standing_filter = np.zeros(z_anom.shape, dtype='complex128')
z_anom_travelling_filter = np.zeros(z_anom.shape, dtype='complex128')

for i in range(z_anom.shape[0]):
    sys.stdout.write('Computing FFT '+str(i+1)+'/'+str(z_anom.shape[0])+'\r')
    sys.stdout.flush()
    z_anom_trans, z_anom_trans_standing, z_anom_trans_travelling = \
            wnfreq.calc_wnfreq_spectrum(z_anom[i], wn_max)
    z_anom_filter[i] = wnfreq.invert_wnfreq_spectrum(
            z_anom_trans*z_filter,1,wn_max,z_lon.size,tol=1e6)
    del z_anom_trans
    gc.collect()

    z_anom_standing[i] = wnfreq.invert_wnfreq_spectrum(
            z_anom_trans_standing,1,wn_max,z_lon.size,tol=1e6)
    z_anom_standing_filter[i] = wnfreq.invert_wnfreq_spectrum(
            z_anom_trans_standing*z_filter,1,wn_max,z_lon.size,tol=1e6)
    del z_anom_trans_standing
    gc.collect()

    z_anom_travelling[i] = wnfreq.invert_wnfreq_spectrum(
            z_anom_trans_travelling,1,wn_max,z_lon.size,tol=1e6)
    z_anom_travelling_filter[i] = wnfreq.invert_wnfreq_spectrum(
            z_anom_trans_travelling*z_filter,1,wn_max,z_lon.size,tol=1e6)
    del z_anom_trans_travelling
    gc.collect()


# plot composite
z_anom_comp = np.zeros((1,comp_length,1,z_lat.size,z_lon.size),dtype='complex128')
z_anom_comp_standing = np.zeros((1,comp_length,1,z_lat.size,z_lon.size),dtype='complex128')
z_anom_comp_travelling = np.zeros((1,comp_length,1,z_lat.size,z_lon.size),dtype='complex128')

z_anom_shift = np.zeros(z_anom_comp.shape,dtype='complex128')
z_anom_shift_standing = np.zeros(z_anom_comp.shape,dtype='complex128')
z_anom_shift_travelling = np.zeros(z_anom_comp.shape,dtype='complex128')

z_anom_comp_filter = np.zeros((1,comp_length,1,z_lat.size,z_lon.size),dtype='complex128')
z_anom_comp_standing_filter = np.zeros((1,comp_length,1,z_lat.size,z_lon.size),dtype='complex128')
z_anom_comp_travelling_filter = np.zeros((1,comp_length,1,z_lat.size,z_lon.size),dtype='complex128')

z_anom_shift_filter = np.zeros(z_anom_comp.shape,dtype='complex128')
z_anom_shift_standing_filter = np.zeros(z_anom_comp.shape,dtype='complex128')
z_anom_shift_travelling_filter = np.zeros(z_anom_comp.shape,dtype='complex128')

t_anom_comp = np.zeros((1,comp_length,1,t_lat.size,t_lon.size))
t_anom_comp_standing = np.zeros((1,comp_length,1,t_lat.size,t_lon.size))
t_anom_comp_travelling = np.zeros((1,comp_length,1,t_lat.size,t_lon.size))

t_anom_shift = np.zeros(t_anom_comp.shape)

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

        # shift and sum geoheight data
        z_anom_shift *= 0
        z_anom_shift_standing *= 0
        z_anom_shift_travelling *= 0

        z_anom_shift[:,:,:,lat_bound[0]:lat_bound[1],lon_bound[0]:lon_bound[1]] = \
            z_anom[year,day0-comp_length//2:day0+comp_length//2,:,
                    lat_bound[0]:lat_bound[1],lon_bound[0]:lon_bound[1]]
        z_anom_shift_standing[:,:,:,lat_bound[0]:lat_bound[1],lon_bound[0]:lon_bound[1]] = \
            z_anom_standing[year,day0-comp_length//2:day0+comp_length//2,:,
                    lat_bound[0]:lat_bound[1],lon_bound[0]:lon_bound[1]]
        z_anom_shift_travelling[:,:,:,lat_bound[0]:lat_bound[1],lon_bound[0]:lon_bound[1]] = \
            z_anom_travelling[year,day0-comp_length//2:day0+comp_length//2,:,
                    lat_bound[0]:lat_bound[1],lon_bound[0]:lon_bound[1]]
        
        z_anom_comp += z_anom_shift
        z_anom_comp_standing += z_anom_shift_standing
        z_anom_comp_travelling += z_anom_shift_travelling

        # shift and sum filtered geoheight data
        z_anom_shift_filter *= 0
        z_anom_shift_standing_filter *= 0
        z_anom_shift_travelling_filter *= 0

        z_anom_shift_filter[:,:,:,lat_bound[0]:lat_bound[1],lon_bound[0]:lon_bound[1]] = \
            z_anom_filter[year,day0-comp_length//2:day0+comp_length//2,:,
                    lat_bound[0]:lat_bound[1],lon_bound[0]:lon_bound[1]]
        z_anom_shift_standing_filter[:,:,:,lat_bound[0]:lat_bound[1],lon_bound[0]:lon_bound[1]] = \
            z_anom_standing_filter[year,day0-comp_length//2:day0+comp_length//2,:,
                    lat_bound[0]:lat_bound[1],lon_bound[0]:lon_bound[1]]
        z_anom_shift_travelling_filter[:,:,:,lat_bound[0]:lat_bound[1],lon_bound[0]:lon_bound[1]] = \
            z_anom_travelling_filter[year,day0-comp_length//2:day0+comp_length//2,:,
                    lat_bound[0]:lat_bound[1],lon_bound[0]:lon_bound[1]]

        z_anom_comp_filter += z_anom_shift
        z_anom_comp_standing_filter += z_anom_shift_standing
        z_anom_comp_travelling_filter += z_anom_shift_travelling

        # shift and sum temperature data
        t_anom_shift *= 0

        t_anom_shift[:,:,:,t_lat_bound[0]:t_lat_bound[1],t_lon_bound[0]:t_lon_bound[1]] = \
            t_anom[year,day0-comp_length//2:day0+comp_length//2,:,
                    t_lat_bound[0]:t_lat_bound[1],t_lon_bound[0]:t_lon_bound[1]]

        t_anom_comp += t_anom_shift

z_anom_comp /= N_hw
z_anom_comp_standing /= N_hw
z_anom_comp_travelling /= N_hw

z_anom_comp_filter /= N_hw
z_anom_comp_standing_filter /= N_hw
z_anom_comp_travelling_filter /= N_hw

t_anom_comp /= N_hw
t_anom_comp_plot = np.zeros((1,)+t_anom_comp.shape)
t_anom_comp_plot[0] = t_anom_comp

del z_anom
del z_anom_standing
del z_anom_travelling

del z_anom_filter
del z_anom_standing_filter
del z_anom_travelling_filter

del z_anom_shift
del z_anom_shift_standing
del z_anom_shift_travelling
gc.collect()
        
"""Plot Composite Geoheight Anomaly Maps"""
plot_data = np.zeros((3,) +\
              z_anom_comp[:, :, :, :, plot_lon_bound[0]:plot_lon_bound[1]].shape, dtype='complex128')
plot_data[0][0] = z_anom_comp[0, :, :, :, plot_lon_bound[0]:plot_lon_bound[1]]
plot_data[1][0] = z_anom_comp_standing[0, :, :, :, plot_lon_bound[0]:plot_lon_bound[1]]
plot_data[2][0] = z_anom_comp_travelling[0, :, :, :, plot_lon_bound[0]:plot_lon_bound[1]]

plot_center = np.tile(center,(3,1))

fig, ax, cax, cbar = hwt.plot_evo(plot_data[:1], 0, comp_length//2, z_lat,
                                z_lon[plot_lon_bound[0]:plot_lon_bound[1]],
                                center=plot_center, days=[0], show_map=False,
                                show_day=False, clev_bound=[-90,90],
                                clev_num=37)
fig.set_size_inches(17,4.5)
ax.set_title('Composite Circulation Across All Heat Wave Start Dates', fontsize=22)
fig.text(0.45, 0.06, 'Longitude', fontsize=18)
fig.text(0.05, 0.5, 'Latitude', rotation='vertical', fontsize=18)
cbar.set_label(r'Geopotential Height Anomaly [m]', fontsize=14)

fig.savefig('figures/composite_geoheight_map.png', format='png')

fig, ax, cax, cbar = hwt.plot_evo(plot_data[:1], 0, comp_length//2, z_lat,
                                    z_lon[plot_lon_bound[0]:plot_lon_bound[1]],
                                    center=plot_center, days=[-6,-3,0,3,6],
                                    show_map=False, clev_bound=[-90,90],
                                    clev_num=37)
fig.set_size_inches(13,13)
ax[0].set_title('Composite Geopotential Height Anomaly Evolution', fontsize=22)
cbar.set_label(r'Geopotential Height Anomaly [m]', fontsize=18)

fig.savefig('figures/composite_geoheight_maps.png', format='png')

fig, ax, cax, cbar = hwt.plot_evo(plot_data, 0, comp_length//2, z_lat,
                                    z_lon[plot_lon_bound[0]:plot_lon_bound[1]],
                                    center=plot_center, days=[-6,-3,0,3,6],
                                    show_map=False, clev_bound=[-90,90],
                                    clev_num=37)
fig.set_size_inches(17, 8.5)
ax[0,0].set_title('Total Geopotential Height Anomaly', fontsize=16)
ax[0,1].set_title('Standing Geopotential Height Anomaly', fontsize=16)
ax[0,2].set_title('Travelling Geopotential Height Anomaly', fontsize=16)
cbar.set_label(r'Geopotential Height Anomaly [m]', fontsize=16)

fig.savefig('figures/composite_geoheight_stand_trav_maps.png', format='png')

""""""

"""Plot Composite Meridional Mean Maps"""
t_plot_data = np.zeros((1,) + t_anom_comp[:, :, :, :, :].shape)
t_plot_data[0][0] = t_anom_comp[0,:,:,:,:]

fig, ax1, cax1, cbar1, ax2, cax2, cbar2 = hwt.plot_custom_meridional_mean_evo(
                plot_data[:1], z_lon[plot_lon_bound[0]:plot_lon_bound[1]], 0,
                comp_length//2, data2=t_plot_data, lon2=t_lon,
                center_lon=plot_center[:,1],
                day_range=[-comp_length//2,comp_length//2],
                clev_bound1=[-90,90], clev_num1=37, clev_bound2=[-2,2],
                clev_num2=21, cax_gap=0.04)

fig.set_size_inches(11,17)
fig.suptitle('Meridional Mean Evolution')
ax1[0].set_title('Composite Geopotential Height Anomaly', fontsize=16)
ax2[0].set_title('Composite Temperature Anomaly', fontsize=16)
fig.text(0.4, 0.055, 'Longitude', fontsize=16)
ax1[0].set_ylabel('Days +/- Heat Wave Start', fontsize=16)
cbar1.set_label(r'Meridional Mean Geopotential Height Anomaly [m]', fontsize=16)
cbar2.set_label(r'Temperature Anomaly [$^\circ$C]', fontsize=16)

fig.savefig('figures/composite_geoheight_temp_evolution_mean.png', format='png')
plt.show()
""""""
