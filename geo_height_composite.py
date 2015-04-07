"""
File: geo_height_composite.py
Author: Charlie White
Email: charlie.white@mail.utoronto.ca
Github: echarliewhite
Description: This script plots the evolution of the standing and travelling
components of geopotential height for the composite of all identified heatwaves
in the given data set.
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
    'Z_GDS0_ISBL', pbound=(300.0,None), latbound=(70.5,10.5))

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
z_anom_standing = np.zeros(z_anom.shape, dtype='complex128')
z_anom_travelling = np.zeros(z_anom.shape, dtype='complex128')

for i in range(z_anom.shape[0]):
    sys.stdout.write('Computing FFT '+str(i+1)+'/'+str(z_anom.shape[0])+'...   ')
    sys.stdout.flush()
    z_anom_trans, z_anom_trans_standing, z_anom_trans_travelling = \
            wnfreq.calc_wnfreq_spectrum(z_anom[i], wn_max)
    sys.stdout.write('Done\n')
    del z_anom_trans
    gc.collect()

    sys.stdout.write('Computing Inverse FFT of Standing Component '+str(i+1)+'/'+str(z_anom.shape[0])+'...   ')
    sys.stdout.flush()
    z_anom_standing[i] = wnfreq.invert_wnfreq_spectrum(
            z_anom_trans_standing,1,wn_max,z_lon.size,tol=1e6)
    sys.stdout.write('Done\n')
    del z_anom_trans_standing
    gc.collect()

    sys.stdout.write('Computing Inverse FFT of Travelling Component '+str(i+1)+'/'+str(z_anom.shape[0])+'...   ')
    sys.stdout.flush()
    z_anom_travelling[i] = wnfreq.invert_wnfreq_spectrum(
            z_anom_trans_travelling,1,wn_max,z_lon.size,tol=1e6)
    sys.stdout.write('Done\n')
    del z_anom_trans_travelling
    gc.collect()

# create composite geopotential for heat waves
comp_length = 30 # number of days centered on day 0 (even # only)

z_anom_comp = np.zeros((1,comp_length,1,z_lat.size,z_lon.size),
                        dtype='complex128')
z_anom_comp_standing = np.zeros((1,comp_length,1,z_lat.size,z_lon.size),
                        dtype='complex128')
z_anom_comp_travelling = np.zeros((1,comp_length,1,z_lat.size,z_lon.size),
                        dtype='complex128')
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
center_ind[0] = (np.searchsorted(z_lat,center[0],sorter=z_lat_sort) - z_lat.size)*-1
center_ind[1] = np.searchsorted(z_lon,center[1],sorter=z_lon_sort)

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

        z_anom_comp += z_anom_shift
        z_anom_comp_standing += z_anom_shift_standing
        z_anom_comp_travelling += z_anom_shift_travelling

del z_anom
del z_anom_standing
del z_anom_travelling
del z_anom_shift
del z_anom_shift_standing
del z_anom_shift_travelling
gc.collect()
        

z_anom_comp /= N_hw
z_anom_comp_standing /= N_hw
z_anom_comp_travelling /= N_hw

# store all 3 data sets in one array for plotting function
# from west pacific to east atlantic
# west_lim = np.where(z_lon==-180)[0][0]
# east_lim = np.where(z_lon==-30.)[0][0]
# south_lim = np.where(z_lat==10.5)[0][0]
# north_lim = np.where(z_lat==70.5)[0][0]
z_anom_comp_all = np.zeros((3,)+z_anom_comp.shape,dtype='complex128')
# z_anom_comp_all = z_anom_comp_all[:,:,:,:,north_lim:south_lim,west_lim:east_lim]
# z_anom_comp_all[0] = z_anom_comp[:,:,:,north_lim:south_lim,west_lim:east_lim]
# z_anom_comp_all[1] = z_anom_comp_standing[:,:,:,north_lim:south_lim,west_lim:east_lim]
# z_anom_comp_all[2] = z_anom_comp_travelling[:,:,:,north_lim:south_lim,west_lim:east_lim]

plot_lon_bound = [np.searchsorted(z_lon,-180,sorter=z_lon_sort),
                 np.searchsorted(z_lon,-30,sorter=z_lon_sort)]

z_anom_comp_all = z_anom_comp_all[:,:,:,:,inner_lat_bound[0]:inner_lat_bound[1],
                                  plot_lon_bound[0]:plot_lon_bound[1]]
z_anom_comp_all[0] = z_anom_comp[:,:,:,inner_lat_bound[0]:inner_lat_bound[1],
                                 plot_lon_bound[0]:plot_lon_bound[1]]
z_anom_comp_all[1] = z_anom_comp_standing[:,:,:,inner_lat_bound[0]:inner_lat_bound[1],
                                          plot_lon_bound[0]:plot_lon_bound[1]]
z_anom_comp_all[2] = z_anom_comp_travelling[:,:,:,inner_lat_bound[0]:inner_lat_bound[1],
                                            plot_lon_bound[0]:plot_lon_bound[1]]

# create evolution plots for composite heat wave
days = [-6,-3,0,3,6]

# fig, ax = hwt.plot_evo(z_anom_comp_all,0,comp_length/2,\
#         z_lat[north_lim:south_lim],z_lon[west_lim:east_lim],
#         center,days,show_map=False)
fig, ax = hwt.plot_evo(z_anom_comp_all,0,comp_length/2,\
        z_lat[inner_lat_bound[0]:inner_lat_bound[1]],z_lon[plot_lon_bound[0]:plot_lon_bound[1]],
        center,days,show_map=False)

ax[0,0].set_title('Composite Total Geoheight')
ax[0,1].set_title('Composite Standing Geoheight')
ax[0,2].set_title('Composite Travelling Geoheight')

#fig.savefig('output/evo_composite.png',format='png')

plt.show()
