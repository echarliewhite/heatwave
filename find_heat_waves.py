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


def import_temps(directory):
    """This function will import the temperature data in the given directory,
    put it all in one array, and return it along with arrays containing
    longitude, latitude, pressure level, and date/time information."""
    
    files = os.listdir(directory)
    var_name='T_GDS0_ISBL' #Temperature

    # determine the shape of the data
    nc = nc4.Dataset(directory+files[0])
    data = nc.variables[var_name][:]
    time = nc.variables['initial_time0_hours'][:]
    plev = nc.variables['lv_ISBL1'][:]
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
        data[i] = nc.variables[var_name][:]
        nc.close()

    # add in year as a returned variable
    return data, time, plev, lat, lon
    
def find_heat_waves(temps,threshold,year):
    """
    Finds heat waves in a given year based on a historical threshold
    percentile at each date and location. Heat wave day is defined as one in
    which the threshold temp is exceeded on that day and the following four
    days in at least 10 locations/grid points. Assumes that in temps the first
    axis [0] represents the year and the second axis [1] represents the date.
    """
    
    # find threshold percentile value at each location/date
    # (year axis removed)
    temps_threshold = np.percentile(temps,threshold,axis=0)
    hot_day_loc = np.greater(temps[year], temps_threshold)
    hot_day = np.empty(temps.shape[1],dtype=bool)
    for i in range(hot_day.size):
        Nhot = np.flatnonzero(hot_day_loc[i]).size
        if Nhot >= 10:
            hot_day[i] = True
        else:
            hot_day[i] = False

    heatwave = np.empty(hot_day.size,dtype=bool)
    heatwave *= 0
    for i in range(heatwave.size-4):
        if hot_day[i]:
            heatwave[i] = True
            for j in range(4):
                if not hot_day[i+j+1]:
                    heatwave[i] = False

    return heatwave
