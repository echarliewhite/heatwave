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
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

z_data, z_time, z_plev, z_lat, z_lon = hwt.import_nc_dir('/Users/charliewhite/Documents/Year4/Thesis/climdata/ERAInterim/dailymean/z/midlatitude/', 'Z_GDS0_ISBL', 300.0)

t_data, t_time, t_plev, t_lat, t_lon = hwt.import_nc_dir('/Users/charliewhite/Documents/Year4/Thesis/climdata/ERAInterim/dailymean/t/summersurface/', 'T_GDS0_ISBL', 1000.0)

# get heat wave data
N_points = 0.05*t_lat.size*t_lon.size # 5% of grid points
heat_wave_dict = hwt.find_heat_waves(t_data,97.5,t_lat,t_lon,N_points,5)

z_clim = np.average(z_data, axis=0) # climatology
# geoheight anomaly
z_anom = z_data
for i in range(z_data.shape[0]):
    z_anom[i] -= z_clim

wn_max = 15
z_anom_trans, z_anom_trans_standing, z_anom_trans_travelling = \
        wnfreq.calc_wnfreq_spectrum(z_anom, wn_max)
z_anom_standing = wnfreq.invert_wnfreq_spectrum(
        z_anom_trans_standing,1,wn_max,z_lon.size,tol=1e6)
z_anom_travelling = wnfreq.invert_wnfreq_spectrum(
        z_anom_trans_travelling,1,wn_max,z_lon.size,tol=1e6)

# day0 for a heat wave
hwave = heat_wave_dict.values()[1]
min_lat_ind = np.where(z_lat==np.min(t_lat))[0][0]
max_lat_ind = np.where(z_lat==np.max(t_lat))[0][0]
min_lon_ind = np.where(z_lon==np.min(t_lon))[0][0]
max_lon_ind = np.where(z_lon==np.max(t_lon))[0][0]
z_anom_standing_day0 = z_anom_standing[hwave.year,hwave.start,0,
        :,min_lon_ind:max_lon_ind+1]
z_anom_travelling_day0 = z_anom_travelling[hwave.year,hwave.start,0,
        :,min_lon_ind:max_lon_ind+1]

# llcrnrlat,llcrnrlon,urcrnrlat,urcrnrlon
# are the lat/lon values of the lower left and upper right corners
# of the map.
# resolution = 'c' means use crude resolution coastlines.
m = Basemap(projection='cyl',llcrnrlat=np.min(t_lat),urcrnrlat=np.max(t_lat),\
            llcrnrlon=np.min(t_lon),urcrnrlon=np.max(t_lon),resolution='c')
# contour levels
clevs = np.arange(np.min(z_anom_standing),np.max(z_anom_standing),
        (np.max(z_anom_standing) - np.min(z_anom_standing))/100.0)
# x,y grid from lat, lon
lons, lats = np.meshgrid(t_lon,t_lat)
x, y = m(lons, lats)
# figure
fig, ax = plt.subplots(1,2)
# contours
plt.sca(ax[0])
cs = m.contour(x,y,z_anom_standing_day0,clevs,colors='k',linewidths=1.)
plt.clabel(cs, inline=1, fontsize=10)
m.drawcoastlines(linewidth=1.25)
m.fillcontinents(color='0.5')
# draw parallels and meridians.
m.drawparallels(np.arange(np.min(t_lat),np.max(t_lat),5.),labels=[1,1,0,0,])
m.drawmeridians(np.arange(np.min(t_lon),np.max(t_lon),10.),labels=[0,0,0,1])
m.drawmapboundary()
plt.title("Day 0 Standing Geo Height Anomaly")

plt.sca(ax[1])
cs = m.contour(x,y,z_anom_travelling_day0,clevs,colors='k',linewidths=1.)
plt.clabel(cs, inline=1, fontsize=10)
m.drawcoastlines(linewidth=1.25)
m.fillcontinents(color='0.5')
# draw parallels and meridians.
m.drawparallels(np.arange(np.min(t_lat),np.max(t_lat),5.),labels=[1,1,0,0,])
m.drawmeridians(np.arange(np.min(t_lon),np.max(t_lon),10.),labels=[0,0,0,1])
m.drawmapboundary()
plt.title("Day 0 Travelling Geo Height Anomaly")
plt.show()

import pdb; pdb.set_trace()
