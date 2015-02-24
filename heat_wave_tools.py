"""
File: find_heat_waves.py
Author: Charlie White
Email: charlie.white@mail.utoronto.ca
Github: echarliewhite
Description: This script contains functions that can identify heat waves based
on a threshold temperature defined relative to historical temperatures.
"""
import numpy as np
import netCDF4 as nc4
import os
import string


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
        

def import_nc_dir(directory, var_name, pressure):
    """
    This function will import the data in the netcdf files in the given
    directory, put it all in one array, and return it along with arrays
    containing longitude, latitude, pressure level, and date/time information.
    The variable of interest (e.g. geo height, temperature) must be specified
    according to its name in the .nc file.
    """
    
    files = os.listdir(directory)
    # don't try to import non-netCDF files
    for item in files:
        if string.find(item, '.nc') == -1:
            files.remove(item)

    # determine the shape of the data
    nc = nc4.Dataset(directory+files[0])
    plev = nc.variables['lv_ISBL1'][:]
    p_ind = np.where(plev==pressure)[0]
    data = nc.variables[var_name][:][:,p_ind,:,:]
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
        data[i] = nc.variables[var_name][:][:,p_ind,:,:]
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

