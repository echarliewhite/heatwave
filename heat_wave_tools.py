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

def plot_evo(data, year, day0, lat, lon, center=None, days=[0], show_map=True,
             cmap=plt.cm.RdBu_r,clev_bound=None, clev_num=20, show_day=True,
             lat_ticks=10, lon_ticks=30):
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
    fig, ax = plt.subplots(len(days), data.shape[0],
                            subplot_kw={'frame_on':True}, facecolor='w',
                            figsize=(17,11))
    marker_style = dict(marker='o', color='black', markersize=10)

    # create basemap instance
    m = Basemap(projection='cyl',llcrnrlat=np.min(lat),\
                urcrnrlat=np.max(lat),llcrnrlon=np.min(lon),\
                urcrnrlon=np.max(lon),resolution='c')
    # make x,y grid from lat,lon
    lon_grid,lat_grid = np.meshgrid(lon,lat)
    x,y = m(lon_grid,lat_grid)
    # contour levels
    if clev_bound is None:
        # extract data in the range of interest for generating contours
        if len(days) == 1:
            data_chunk = np.real(data[:,year,day0+days[0]])
        else:
            data_chunk = np.real(data[:,year,day0+days[0]:day0+days[-1]])
        extreme = np.max(np.absolute(data_chunk))
        clevs = np.linspace(-((extreme+10)//10)*10,((extreme+10)//10)*10,
                            num=clev_num)
    else:
        clevs = np.linspace(clev_bound[0],clev_bound[1],num=clev_num)

    # fill in individual plots
    for i in range(data.shape[0]):
        for j in range(len(days)):
            # set current axes
            if len(days) == 1:
                if data.shape[0] == 1:
                    plt.sca(ax)
                else:
                    plt.sca(ax[i])
            else:
                if data.shape[0] == 1:
                    plt.sca(ax[j])
                else:
                    plt.sca(ax[j,i])

            # check that the day is in range
            if 0 <= day0 + days[j] < data.shape[2]:
                plot_data = np.real(data[i][year][day0+days[j]][0])
                if show_map:
                    # draw map features
                    m.drawcoastlines(linewidth=2)
                    m.drawcountries(linewidth=2)

                if i == 0:
                    m.drawparallels(np.arange((np.min(lat)//10)*10+10,np.max(lat),lat_ticks),labels=[1,0,0,0])
                else:
                    m.drawparallels(np.arange((np.min(lat)//10)*10+10,np.max(lat),lat_ticks))

                if j == len(days) - 1:
                    m.drawmeridians(np.arange((np.min(lon)//10)*10+10,np.max(lon),lon_ticks),labels=[0,0,0,1])
                else:
                    m.drawmeridians(np.arange((np.min(lon)//10)*10+10,np.max(lon),lon_ticks))
                # fit plot to axis
                # plot contours
                cs = m.contourf(x,y,plot_data,clevs,cmap=cmap)
                csf = m.contour(x,y,plot_data,clevs,linewidths=1,colors='k')
                # label day number
                if show_day == True:
                    day_lab = AnchoredText('Day %d' % days[j], loc=1, frameon=True)
                    day_lab.patch.set_boxstyle('round,pad=0.,rounding_size=0.2')
                    if data.shape[0] == 1:
                        ax[j].add_artist(day_lab)
                    else:
                        ax[j,i].add_artist(day_lab)
                # plot center point
                if center is not None:
                    plt.plot(center[i,1],center[i,0], **marker_style)
            else:
                plt.text(0.2,0.5,'Date out of range')
    plt.subplots_adjust(left=0.1, right=0.85, wspace=0.01, hspace=0.05)
    cax = fig.add_axes([0.875, 0.1, 0.025, 0.8])
    cbar = plt.colorbar(cs, cax=cax, orientation='vertical')
    return fig, ax, cax, cbar

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

    z_cs1 = z_ax1.contourf(zx,zy,z_anom_mer_mean,z_clevs,cmap=plt.cm.RdBu_r)
    t_cs = t_ax.contourf(tx,ty,t_anom_mer_mean,t_clevs,cmap='rainbow')
    z_cs2 = z_ax2.contourf(zx,zy,z_anom_standing_mer_mean,z_clevs,cmap=plt.cm.RdBu_r)
    z_cs3 = z_ax3.contourf(zx,zy,z_anom_travelling_mer_mean,z_clevs,cmap=plt.cm.RdBu_r)

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

def plot_custom_meridional_mean_evo(data1,  lon1, year,
        day0, data2=None, lon2=None, center_lon=None, day_range=[-30,30],
        cmap1=plt.cm.RdBu_r, cmap2='rainbow', clev_bound1=None, clev_num1=20,
        clev_bound2=None, clev_num2=20, cax_gap=0.025):
    """
    This function will create a set of plots showing the day by day evolution
    of the meridional average of the given data.Fot reference see fig. 1 of
    Kushner/Watt-Meyer's decomposition paper.

    data - array containing z data with dimensions [# of
             datasets,years,days,plev=1,z_lat,z_lon]
    lon - coordinate information for corresponding array
    year - the year in which the dates to be plotted reside
    day0 - heat wave start date, or other date upon which to center the plots
    center_lon - a longitude to identify with a line on the plots
    day_range - number of days to plot on either side of day0.

    Returns:
    fig - a pyplot figure
    ax - an array pyplot axes
    cax - colorbar axes
    """
    # create figure/axes
    fig = plt.figure(figsize=(17,11), facecolor='w')

    left = 0.1
    bottom = 0.1
    width = 0.8
    height = 0.8
    gap = 0.01
    cax_width = 0.025

    if data2 is None:
        total_lon = lon1.size*data1.shape[0]
        total_plot_width = width - (cax_width+cax_gap) - gap*data1.shape[0]
        panel_width1 = (lon1.size/float(total_lon))*total_plot_width
    else:
        total_lon = lon1.size*data1.shape[0] + lon2.size*data2.shape[0]
        total_plot_width = width - (cax_width+cax_gap)*2 - gap*(data1.shape[0]+data2.shape[0])
        panel_width1 = (lon1.size/float(total_lon))*total_plot_width
        panel_width2 = (lon2.size/float(total_lon))*total_plot_width

    # axis boundaries
    rect1 = np.zeros((data1.shape[0],4))
    ax1 = np.empty(data1.shape[0], dtype=object)
    for i in range(data1.shape[0]):
        rect1[i][0] = left + (panel_width1+gap)*i
        rect1[i][1] = bottom
        rect1[i][2] = panel_width1
        rect1[i][3] = height
        
        if i == 0:
            ax1[i] = fig.add_axes(rect1[i])
        else:
            ax1[i] = fig.add_axes(rect1[i],sharey=ax1[0])
    
    # colorbar axes
    cax_rect1 = np.zeros(4)
    cax_rect1[0] = left + (panel_width1+gap)*data1.shape[0]
    cax_rect1[1] = bottom
    cax_rect1[2] = cax_width
    cax_rect1[3] = height

    cax1 = fig.add_axes(cax_rect1)

    if data2 is not None:
        # axis boundaries
        rect2 = np.zeros((data2.shape[0],4))
        ax2 = np.empty(data2.shape[0], dtype=object)
        for i in range(data2.shape[0]):
            rect2[i][0] = left + (panel_width1+gap)*data1.shape[0] +\
                (cax_width+cax_gap) + (panel_width2+gap)*i
            rect2[i][1] = bottom
            rect2[i][2] = panel_width2
            rect2[i][3] = height
            
            ax2[i] = fig.add_axes(rect2[i],sharey=ax1[0])
    
        # colorbar axes
        cax_rect2 = np.zeros(4)
        cax_rect2[0] = left + (panel_width1+gap)*data1.shape[0] +\
            (cax_width+cax_gap) + (panel_width2+gap)*data2.shape[0]
        cax_rect2[1] = bottom
        cax_rect2[2] = cax_width
        cax_rect2[3] = height

        cax2 = fig.add_axes(cax_rect2)

    # take meridional averages and slice out time range of interest
    data1_mer_mean = np.average(np.real(data1),axis=4)[:,year,
                                    day0+day_range[0]:day0+day_range[1],:,:]
    if data2 is not None:
        data2_mer_mean = np.average(np.real(data2),axis=4)[:,year,
                                    day0+day_range[0]:day0+day_range[1],:,:]

    # set up contour plot
    days = np.arange(day_range[0],day_range[1])
    x1,y1 = np.meshgrid(lon1,days)
    if data2 is not None:
        x2,y2 = np.meshgrid(lon2,days)

    if clev_bound1 is None:
        extreme = np.max(np.absolute(data1_mer_mean))
        clevs1 = np.linspace(-((extreme+10)//10)*10,
                            ((extreme+10)//10)*10,num=clev_num1)
    else:
        clevs1 = np.linspace(clev_bound1[0],clev_bound1[1],num=clev_num1)

    if data2 is not None:
        if clev_bound2 is None:
            extreme = np.max(np.absolute(data2_mer_mean))
            clevs2 = np.linspace(-((extreme+10)//10)*10,
                            ((extreme+10)//10)*10,num=clev_num2)
        else:
            clevs2 = np.linspace(clev_bound2[0],clev_bound2[1],num=clev_num2)
    cs1 = np.empty(data1.shape[0],dtype=object)
    for i in range(data1.shape[0]):
        cs1[i] = ax1[i].contourf(x1,y1,np.squeeze(data1_mer_mean[i]),clevs1,cmap=cmap1)
        cs = ax1[i].contour(x1,y1,np.squeeze(data1_mer_mean[i]), clevs1,
                                linewidths=1, colors='k')
    cbar1 = plt.colorbar(cs1[0], cax=cax1, orientation='vertical')

    if data2 is not None:
        cs2 = np.empty(data2.shape[0],dtype=object)
        for i in range(data2.shape[0]):
            cs2[i] = ax2[i].contourf(x2,y2,np.squeeze(data2_mer_mean[i]),clevs2,cmap=cmap2)
        cbar2 = plt.colorbar(cs2[0], cax=cax2, orientation='vertical')

    # identify center and day0
    # format ticks
    ax1[0].yaxis.set_major_locator(plt.MultipleLocator(3.0))
    for i in range(data1.shape[0]):
        ax1[i].xaxis.set_major_locator(plt.MultipleLocator(30.0))
        if i > 0:
            plt.setp(ax1[i].get_yticklabels(), visible=False)
        ax1[i].hlines(0,lon1[0],lon1[-1],linestyles='dashed')
        if center_lon is not None:
            ax1[i].vlines(center_lon[i],days[0],days[-1],linestyles='dashed')

    if data2 is not None:
        for i in range(data2.shape[0]):
            ax2[i].xaxis.set_major_locator(plt.MultipleLocator(30.0))
            plt.setp(ax2[i].get_yticklabels(), visible=False)
            ax2[i].hlines(0,lon2[0],lon2[-1],linestyles='dashed')
            if center_lon is not None:
                ax2[i].vlines(center_lon[0],days[0],days[-1],linestyles='dashed')

    if data2 is not None:
        return fig, ax1, cax1, cbar1, ax2, cax2, cbar2
    else:
        return fig, ax1, cax1, cbar1
