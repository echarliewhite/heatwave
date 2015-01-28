import find_heat_waves as fhw
import matplotlib.pyplot as plt
import matplotlib.collections as col
import numpy as np

# arguments
threshold = 97.5
directory = '/Users/charliewhite/Documents/Year4/Thesis/climdata/ERAInterim/dailymean/t/summersurface/'

# find heat waves in given year
temps, time, plev, lat, lon = fhw.import_temps(directory)
N_hw_days = np.zeros(temps.shape[0])
for j in range(temps.shape[0]):
    heatwavedays = fhw.find_heat_waves(temps,threshold,j,lat,lon,44,5)
    N_hw_days[j] = np.count_nonzero(heatwavedays)

    # what is the daily average temperature in the given year, and across all years
    avtemp = np.zeros(time.size)
    for i in range(time.size):
        avtemp[i] = np.average(temps[:,i])
    yearavtemp = np.zeros(time.size)
    for i in range(time.size):
        yearavtemp[i] = np.average(temps[j][i])

    # plot
    if j == 1995 - 1979:
        fig, ax = plt.subplots(1,1)
        ax.plot(time,yearavtemp,'b')
        ax.plot(time,avtemp,'g')
        collection = col.BrokenBarHCollection.span_where(time,where=heatwavedays,
                ymin=min(yearavtemp),ymax=max(yearavtemp),facecolor='red',alpha=0.5)
        ax.add_collection(collection)


fig, ax = plt.subplots(1,1)
ax.bar(np.arange(temps.shape[0])+1979,N_hw_days)

print np.sum(N_hw_days)
plt.show()
        
