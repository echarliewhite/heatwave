"""
File: find_heat_waves.py
Author: Charlie White
Email: charlie.white@mail.utoronto.ca
Github: echarliewhite
Description: This module contains functions that can identify heat waves based
on a threshold temperature defined relative to historical temperatures. A class
for storing heat waves and plotting functions are included as well.
"""
import numpy as np
import netCDF4 as nc4
import os
import sys
import string
from mpl_toolkits.basemap import Basemap
from matplotlib.offsetbox import AnchoredOffsetbox, AnchoredText, TextArea
import matplotlib.pyplot as plt

class heat_wave:
    """
    A class for storing heatwaves and their properties.
    """
    def __init__(self, year, start, length, lat, lon):
        self.year = year # year of the heat wave
        self.start = start # starting date of the heat wave
        self.length = length # number of days the heat wave lasted
        self.lat = lat # average latitude of grid points above temp threshold
        self.lon = lon # average longitude of grid points above temp threshold
        
class AnchoredText(AnchoredOffsetbox):
    """
    for adding labels within figures
    """
    def __init__(self, s, loc, pad=0.4, borderpad=0.5, prop=None, frameon=True):
        self.txt = TextArea(s,minimumdescent=False)

        super(AnchoredText, self).__init__(loc, pad=pad, borderpad=borderpad,
                                       child=self.txt,
                                       prop=prop,
                                       frameon=frameon)

def import_nc_dir(directory, var_name, pbound=(None, None),
                    latbound=(None,None), lonbound=(None,None)):
    """
    This function will import the data in the netcdf files in the given
    directory, put it all in one array, and return it along with arrays
    containing longitude, latitude, pressure level, and date/time information.
    The variable of interest (e.g. geo height, temperature) must be specified
    according to its name in the .nc file. To reduce the amount of data
    imported, a single pressure level or a range must be specified.
    """
    
    files = os.listdir(directory)
    # don't try to import non-netCDF files
    for item in files:
        if string.find(item, '.nc') == -1:
            files.remove(item)

    # determine the shape of the data
    nc = nc4.Dataset(directory+files[0])

    # get pressure range
    plev = nc.variables['lv_ISBL1'][:]
    if pbound[0] is None:
        p_ind1 = 0
        p_ind2 = -1
    else:
        p_ind1 = np.where(plev==pbound[0])[0][0]
        if pbound[1] is None:
            p_ind2 = p_ind1 + 1
        else:
            p_ind2 = np.where(p==pbound[1])[0][0] + 1
    plev = plev[p_ind1:p_ind2]

    # have to add hours/year after year 0 to get absolute time
    time = nc.variables['initial_time0_hours'][:] 

    # get latitude range
    lat = nc.variables['g0_lat_2'][:]
    if latbound[0] is None:
        lat_ind1 = 0
        lat_ind2 = -1
    else:
        lat_ind1 = np.where(lat==latbound[0])[0][0]
        if latbound[1] is None:
            lat_ind2 = lat_ind1 + 1
        else:
            lat_ind2 = np.where(lat==latbound[1])[0][0] + 1
    lat = lat[lat_ind1:lat_ind2]

    # get longitude range
    lon = nc.variables['g0_lon_3'][:]
    if lonbound[0] is None:
        lon_ind1 = 0
        lon_ind2 = -1
    else:
        lon_ind1 = np.where(lon==lonbound[0])[0][0]
        if lonbound[1] is None:
            lon_ind2 = lon_ind1 + 1
        else:
            lon_ind2 = np.where(lon==lonbound[1])[0][0] + 1
    lon = lon[lon_ind1:lon_ind2]

    data = nc.variables[var_name][:][:,p_ind1:p_ind2,lat_ind1:lat_ind2,
                                     lon_ind1:lon_ind2]
    # add dimension for years
    Nyears = len(files)
    data = np.zeros((Nyears,)+data.shape)
    # reminder: add variable for years and extract starting year from data
    nc.close()

    # extract data for each year and put in one array
    for i in range(Nyears):
        nc = nc4.Dataset(directory+files[i])
        data[i] = nc.variables[var_name][:][:,p_ind1:p_ind2,
                                            lat_ind1:lat_ind2,lon_ind1:lon_ind2]
        nc.close()
        sys.stdout.write(var_name + ' Import: %d%% Complete\r' % ((i+1)*100/Nyears))
        sys.stdout.flush()
    sys.stdout.write(var_name + ' Import: %d%% Complete\n' % ((i+1)*100/Nyears))
    sys.stdout.flush()

    # add in year as a returned variable
    return data, time, plev, lat, lon
    

def find_heat_waves(temps,threshold,lat,lon,N_points,deg_move):
    """
    Finds heat waves based on a historical threshold percentile at each date
    and location. A heat wave day is defined as one in which the threshold temp
    is exceeded on that day and the following four days in at least N_points
    locations/grid points. Assumes that temps has structure
    [year,day,plev,lat,lon]. lat, lon define the latitudes/longitudes
    represented in temps. N_points determines the number of points which must
    exceed the threshold temperature. deg_move determines how far the center of
    the hot locations may move from day to day in degrees latitude OR longitude
    (to prevent e.g. east coast heat followed coincidentally by west coast heat
    from being considered a single event.)

    Returns a dictionary containing all heat waves. Each key is a heat wave
    class, each entry may be used to contain plots of that heat wave.
    """
    
    # find threshold percentile value at each location/date
    # (year axis removed)
    temps_threshold = np.zeros(temps.shape[1:])
    # look at all temperatures in a 15-day window
    for i in range(temps_threshold.shape[0]):
        # ranges change near ends of array
        lowbound = min(i,7)
        highbound = min(temps_threshold.shape[0] - 1 - i, 7)
        temps_threshold[i,:,:,:] = \
          np.percentile(temps[:,i-lowbound:i+highbound,:,:,:],threshold,axis=(0,1))

    heat_wave_dict = dict()

    # find days that exceed the threshold value
    hot_day_loc = np.empty(temps.shape, dtype=bool)
    for i in range(temps.shape[0]): # iterate through years
        # days and locations exceeding threshold temperature
        hot_day_loc[i] = np.greater(temps[i], temps_threshold)
        prev_day_heat_wave = False # bool indicating whether the day before was "hot"
        for j in range(temps.shape[1]): # iterate through days
            # hot grid points greater than specified threshold?
            if np.flatnonzero(hot_day_loc[i][j]).size >= N_points:
                if not prev_day_heat_wave: # possible heat wave start date
                    # locations where it is hot
                    hot_day_loc_indices = np.nonzero(hot_day_loc[i][j])
                    # find center of hot locations
                    lat_av = np.average(lat[hot_day_loc_indices[1]])
                    lon_av = np.average(lon[hot_day_loc_indices[2]])
                    
                    k = 1
                    if j < temps.shape[1] - 1:
                        next_day_hot = True
                    else: next_day_hot = False

                    while next_day_hot:
                        # is the next day "hot"?
                        if np.flatnonzero(hot_day_loc[i][j+k]).size >= N_points:
                            hot_day_loc_indices = np.nonzero(hot_day_loc[i][j+k])
                            lat_av_next = np.average(lat[hot_day_loc_indices[1]])
                            lon_av_next = np.average(lon[hot_day_loc_indices[2]])
                            # does the center of the heat wave move too quickly?
                            if (abs(lat_av_next - lat_av) > deg_move or
                                  abs(lon_av_next - lon_av) > deg_move):
                                next_day_hot = False
                            else:
                                k += 1
                                if j + k >= temps.shape[1]:
                                    next_day_hot = False
                                lat_av, lon_av = lat_av_next, lon_av_next
                        else: next_day_hot = False

                    if k > 4: # sufficiently long string of hot days
                        hot_day_loc_indices = np.nonzero(hot_day_loc[i][j])
                        lat_av = np.average(lat[hot_day_loc_indices[1]])
                        lon_av = np.average(lon[hot_day_loc_indices[2]])
                        # store heat wave
                        heat_wave_dict['year_%d_day_%d' % (i,j)] = heat_wave(i,j,k,lat_av,lon_av)
                        prev_day_heat_wave = True
                else:
                    k -= 1
                    # when heat wave ends start looking for heat waves again
                    if k == 1: prev_day_heat_wave = False

        sys.stdout.write('Heat Wave Search: %d%% Complete\r' % 
                            ((i+1)*100/temps.shape[0]))
        sys.stdout.flush()
    sys.stdout.write('Heat Wave Search: %d%% Complete\n' % 
                        ((i+1)*100/temps.shape[0]))
    sys.stdout.flush()
    return heat_wave_dict

def plot_evo(data, year, day0, lat, lon, center=None, days=[0], show_map=True):
    """
    This function will create a set of plots (e.g.
    total/standing/travelling) side by side showing the evolution of the data
    sets.

    data - array containing data with dimensions 
            [# of datasets,year,day,plev=1,lat,lon]
    year - year 0 refers to first index in data (e.g. don't use 1982, use 3 if
            data starts in 1979)
    day0 - what day to center the evolution on
    lat - array of latitudes corresponding to data
    lon - array of longitudes corresponding to data
    center - a point to plot on the map (e.g. center of heat wave)
    days - list of integers representing days relative to day0 for which to
            create plots

    Returns:
    fig - a pyplot figure
    ax - an array of pyplot axes with dimensions [days.size,# of datasets]
    """
    # create figure/axes
    fig, ax = plt.subplots(len(days),data.shape[0],figsize=(17,12),
                            sharex=True,sharey=True)

    # create basemap instance
    m = Basemap(projection='cyl',llcrnrlat=np.min(lat),\
                urcrnrlat=np.max(lat),llcrnrlon=np.min(lon),\
                urcrnrlon=np.max(lon),resolution='c')
    # make x,y grid from lat,lon
    lon_grid,lat_grid = np.meshgrid(lon,lat)
    x,y = m(lon_grid,lat_grid)
    # extract data in the range of interest for generating contours
    data_chunk = np.real(data[:,year,day0+days[0]:day0+days[-1]])
    # contour levels
    extreme = np.max(np.absolute(data_chunk))
    clevs = np.linspace(-((extreme+10)//10)*10,((extreme+10)//10)*10,num=20)

    # fill in individual plots
    for i in range(data.shape[0]):
        for j in range(len(days)):
            # set current axes
            plt.sca(ax[j,i])

            # check that the day is in range
            if 0 <= day0 + days[j] < data.shape[2]:
                plot_data = np.real(data[i][year][day0+days[j]][0])
                if show_map:
                    # draw map features
                    m.drawcoastlines(linewidth=1.25)
                m.drawparallels(np.arange(np.min(lat),np.max(lat),5.))
                m.drawmeridians(np.arange(np.min(lon),np.max(lon),10.))
                m.drawmapboundary()
                # fit plot to axis
                # plot contours
                cs = m.contourf(x,y,plot_data,clevs,cmap='bwr')
                cs = m.contour(x,y,plot_data,clevs,linewidths=1,colors='k')
                # label day number
                day_lab = AnchoredText('Day %d' % days[j], loc=1, frameon=True)
                day_lab.patch.set_boxstyle('round,pad=0.,rounding_size=0.2')
                ax[j,i].add_artist(day_lab)
                # plot center point
                if center is not None:
                    plt.plot(center[1],center[0],'g*')
            else:
                plt.text(0.2,0.5,'Date out of range')

    plt.tight_layout() 
    return fig, ax

def plot_meridional_mean_evo(z_data, z_lon, t_data, t_lon, year,
                             day0, center_lon=None, day_range=[-30,30]):
    """
    This function will create a set of plots showing the day by day evolution
    of the meridional average of the total, standing, and travelling geoheight
    anomaly alongside the meridional temperature average. For reference see
    fig. 1 of Kushner/Watt-Meyer's decomposition paper.

    z_data - array containing z data with dimensions [datasets=3,years,days,plev=1,z_lat,z_lon]
             order of datasets - total anomaly, standing, travelling
    t_data - array containing T data with dimensions [years,days,plev=1,t_lat,t_lon]
    z_lon, t_lon - coordinate information for corresponding array
    year - the year in which the dates to be plotted reside
    day0 - heat wave start date, or other date upon which to center the plots
    center_lon - a longitude to identify with a line on the plots
    day_range - number of days to plot on either side of day0.

    Returns:
    fig - a pyplot figure
    ax - an array of 4 pyplot axes
    """
    # create figure/axes
    fig = plt.figure(figsize=(17,11))

    left = 0.1
    bottom = 0.1
    width = 0.8
    height = 0.8

    total_lon_width = z_lon.size*3 + t_lon.size
    z_panel_width = z_lon.size*width/total_lon_width
    t_panel_width = t_lon.size*width/total_lon_width

    z_rect1 = [left, bottom, z_panel_width, height]
    t_rect = [left+z_panel_width, bottom, t_panel_width, height]
    z_rect2 = [left+z_panel_width+t_panel_width, bottom, z_panel_width, height]
    z_rect3 = [left+z_panel_width*2+t_panel_width, bottom, z_panel_width, height]
    
    z_ax1 = fig.add_axes(z_rect1)
    t_ax = fig.add_axes(t_rect, sharey=z_ax1)
    z_ax2 = fig.add_axes(z_rect2, sharey=z_ax1)
    z_ax3 = fig.add_axes(z_rect3, sharey=z_ax1)

    # colorbar axes
    z_cax_rect = [left+width, bottom, (1. - (left+width))/4., height]
    t_cax_rect = [left+width+2*(1. - (left+width))/4., bottom, (1. - (left+width))/4., height]

    z_cax = fig.add_axes(z_cax_rect)
    t_cax = fig.add_axes(t_cax_rect)

    # take meridional averages
    z_data_mer_mean = np.average(np.real(z_data),axis=4)
    t_data_mer_mean = np.average(t_data,axis=3)
    
    z_anom_mer_mean = np.squeeze(z_data_mer_mean[0,year,day0+day_range[0]:day0+day_range[1],:,:])
    t_anom_mer_mean = np.squeeze(t_data_mer_mean[year,day0+day_range[0]:day0+day_range[1],:,:])
    z_anom_standing_mer_mean = np.squeeze(z_data_mer_mean[1,year,day0+day_range[0]:day0+day_range[1],:,:])
    z_anom_travelling_mer_mean = np.squeeze(z_data_mer_mean[2,year,day0+day_range[0]:day0+day_range[1],:,:])

    # set up contour plot
    days = np.arange(day_range[0],day_range[1])
    zx,zy = np.meshgrid(z_lon,days)
    tx,ty = np.meshgrid(t_lon,days)

    z_extreme = np.max(np.absolute(z_data_mer_mean))
    z_clevs = np.linspace(-((z_extreme+10)//10)*10,((z_extreme+10)//10)*10,num=20)
    t_extreme = np.max(np.absolute(t_anom_mer_mean))
    t_clevs = np.linspace(-((t_extreme+2)//2)*2,((t_extreme+2)//2)*2,num=10)

    z_cs1 = z_ax1.contourf(zx,zy,z_anom_mer_mean,z_clevs,cmap='bwr')
    t_cs = t_ax.contourf(tx,ty,t_anom_mer_mean,t_clevs,cmap='autumn')
    z_cs2 = z_ax2.contourf(zx,zy,z_anom_standing_mer_mean,z_clevs,cmap='bwr')
    z_cs3 = z_ax3.contourf(zx,zy,z_anom_travelling_mer_mean,z_clevs,cmap='bwr')

    # identify center and day0
    z_ax1.hlines(0,z_lon[0],z_lon[-1],linestyles='dashed')
    t_ax.hlines(0,t_lon[0],t_lon[-1],linestyles='dashed')
    z_ax2.hlines(0,z_lon[0],z_lon[-1],linestyles='dashed')
    z_ax3.hlines(0,z_lon[0],z_lon[-1],linestyles='dashed')

    if center_lon is not None:
        z_ax1.vlines(center_lon,days[0],days[-1],linestyles='dashed')
        t_ax.vlines(center_lon,days[0],days[-1],linestyles='dashed')
        z_ax2.vlines(center_lon,days[0],days[-1],linestyles='dashed')
        z_ax3.vlines(center_lon,days[0],days[-1],linestyles='dashed')

    # colorbars
    plt.colorbar(z_cs1, cax=z_cax, orientation='vertical')
    plt.colorbar(t_cs, cax=t_cax, orientation='vertical')

    # labels and ticks
    z_ax1.xaxis.set_major_locator(plt.MultipleLocator(30.0))
    t_ax.xaxis.set_major_locator(plt.MultipleLocator(30.0))
    z_ax2.xaxis.set_major_locator(plt.MultipleLocator(30.0))
    z_ax3.xaxis.set_major_locator(plt.MultipleLocator(30.0))

    z_ax1.yaxis.set_major_locator(plt.MultipleLocator(3.0))
    plt.setp(t_ax.get_yticklabels(), visible=False)
    plt.setp(z_ax2.get_yticklabels(), visible=False)
    plt.setp(z_ax3.get_yticklabels(), visible=False)

    z_ax1.set_title('Total Geoheight Anomaly',fontsize=12)
    t_ax.set_title('Temperature Anomaly',fontsize=12)
    z_ax2.set_title('Standing Geoheight Anomaly',fontsize=12)
    z_ax3.set_title('Travelling Geoheight Anomaly',fontsize=12)

    z_ax1.set_xlabel('Longitude')
    t_ax.set_xlabel('Longitude')
    z_ax2.set_xlabel('Longitude')
    z_ax3.set_xlabel('Longitude')

    z_ax1.set_ylabel('Day #')
    
    plt.tick_params(axis='both', which='major', labelsize='8')

    ax = [z_ax1, t_ax, z_ax2, z_ax3]
    return fig, ax
