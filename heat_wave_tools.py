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

def import_nc_dir(directory, var_name, pressure1=1000, pressure2=None):
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
    plev = nc.variables['lv_ISBL1'][:]
    p_ind1 = np.where(plev==pressure1)[0][0] # check that this works
    if pressure2 is None:
        p_ind2 = p_ind1 + 1
    else:
        p_ind2 = np.where(plev==pressure2)[0][0] + 1
    plev = plev[p_ind1:p_ind2]
    data = nc.variables[var_name][:][:,p_ind1:p_ind2,:,:]
    # have to add hours/year after year 0 to get absolute time
    time = nc.variables['initial_time0_hours'][:] 
    lat = nc.variables['g0_lat_2'][:]
    lon = nc.variables['g0_lon_3'][:]

    # add dimension for years
    Nyears = len(files)
    data = np.zeros((Nyears,)+data.shape)
    # add variable for years and extract starting year from data
    nc.close()

    # extract temperature data for each year and put in one array
    for i in range(len(files)):
        nc = nc4.Dataset(directory+files[i])
        data[i] = nc.variables[var_name][:][:,p_ind1:p_ind2,:,:]
        nc.close()

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
    return heat_wave_dict

def plot_evo(data, year, day0, lat, lon, center=None, days=[0]):
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
    fig, ax = plt.subplots(len(days),data.shape[0],figsize=(17,11))

    m = Basemap(projection='cyl',llcrnrlat=np.min(lat),\
                urcrnrlat=np.max(lat),llcrnrlon=np.min(lon),\
                urcrnrlon=np.max(lon),resolution='c')
    # make x,y grid from lat,lon
    lon_grid,lat_grid = np.meshgrid(lon,lat)
    x,y = m(lon_grid,lat_grid)

    # fill in individual plots
    for i in range(data.shape[0]):
        for j in range(len(days)):
            # set current axes
            plt.sca(ax[j,i])

            # check that the day is in range
            if 0 <= day0 + days[j] < data.shape[2]:
                plot_data = np.real(data[i][year][day0+days[j]][0])
                # draw map features
                m.drawcoastlines(linewidth=1.25)
                m.fillcontinents()
                m.drawparallels(np.arange(np.min(lat),np.max(lat),5.))
                m.drawmeridians(np.arange(np.min(lon),np.max(lon),10.))
                m.drawmapboundary()
                # contour levels
                clevs = np.arange(np.min(plot_data),np.max(plot_data),
                            (np.max(plot_data) - np.min(plot_data))/20.0)
                # fit plot to axis
                # plot contours
                cs = m.contour(x,y,plot_data,clevs,linewidths=1,colors='k')
                plt.clabel(cs,inline=1,fontsize=10)
                # label day number
                day_lab = AnchoredText('Day %d' % days[j], loc=1, frameon=True)
                day_lab.patch.set_boxstyle('round,pad=0.,rounding_size=0.2')
                ax[j,i].add_artist(day_lab)
                # plot center point
                if center is not None:
                    plt.plot(center[1],center[0],'ro')
            else:
                plt.text(0.2,0.5,'Date out of range')
    plt.tight_layout() 
    return fig, ax

