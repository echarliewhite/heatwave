"""
File: geo_height_evolution.py
Author: Charlie White
Email: charlie.white@mail.utoronto.ca
Github: echarliewhite
Description: This script plots the evolution of the standing and travelling
components of geopotential height for a specified heat wave.
"""
import numpy as np
import heat_wave_tools as hwt
import wnfreq_routines_2_0 as wnfreq
import matplotlib.pyplot as plt
import sys

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

wn_max = 15
sys.stdout.write('Computing FFT...   ')
sys.stdout.flush()
z_anom_trans, z_anom_trans_standing, z_anom_trans_travelling = \
        wnfreq.calc_wnfreq_spectrum(z_anom, wn_max)
sys.stdout.write('Done\n')
sys.stdout.write('Computing Inverse FFT of Standing Component...   ')
sys.stdout.flush()
z_anom_standing = wnfreq.invert_wnfreq_spectrum(
        z_anom_trans_standing,1,wn_max,z_lon.size,tol=1e6)
sys.stdout.write('Done\n')
sys.stdout.write('Computing Inverse FFT of Travelling Component...   ')
sys.stdout.flush()
z_anom_travelling = wnfreq.invert_wnfreq_spectrum(
        z_anom_trans_travelling,1,wn_max,z_lon.size,tol=1e6)
sys.stdout.write('Done\n')

# store all 3 data sets in one array for plotting function
# from west pacific to east atlantic
west_lim = np.where(z_lon==-180)[0][0]
east_lim = np.where(z_lon==-30.)[0][0] + 1
south_lim = np.where(z_lat==10.5)[0][0] + 1
north_lim = np.where(z_lat==70.5)[0][0]
z_anom_all = np.zeros((3,)+z_anom.shape,dtype='complex128')
z_anom_all = z_anom_all[:,:,:,:,north_lim:south_lim,west_lim:east_lim]
z_anom_all[0] = z_anom[:,:,:,north_lim:south_lim,west_lim:east_lim]
z_anom_all[1] = z_anom_standing[:,:,:,north_lim:south_lim,west_lim:east_lim]
z_anom_all[2] = z_anom_travelling[:,:,:,north_lim:south_lim,west_lim:east_lim]

# create evolution plots for each heat wave
days = [-6,-3,0,3,6]
evo_plot_dict = dict()

# for key in heat_wave_dict.keys():
#     hwave = heat_wave_dict[key]
#     year = hwave.year
#     day0 = hwave.start
#     center = [hwave.lat,hwave.lon]
#     fig, ax = hwt.plot_evo(z_anom_all,year,day0,\
#                             z_lat[west_lim:east_lim],z_lon,center,days)
#     evo_plot_dict[key] = (fig, ax)
#     ax[0,0].set_title(key + 'Total Geoheight')
#     ax[0,1].set_title(key + 'Standing Geoheight')
#     ax[0,2].set_title(key + 'Travelling Geoheight')

key = heat_wave_dict.keys()[0]
hwave = heat_wave_dict[key]
year = hwave.year
day0 = hwave.start
center = [hwave.lat,hwave.lon]
fig, ax = hwt.plot_evo(z_anom_all,year,day0,\
        z_lat[north_lim:south_lim],z_lon[west_lim:east_lim],center,days)
evo_plot_dict[key] = (fig, ax)

ax[0,0].set_title(key + ' Total Geoheight')
ax[0,1].set_title(key + ' Standing Geoheight')
ax[0,2].set_title(key + ' Travelling Geoheight')

plt.show()
