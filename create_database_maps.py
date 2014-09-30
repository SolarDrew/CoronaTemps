# -*- coding: utf-8 -*-
"""
Created on Wed Sep  3 16:13:05 2014

@author: drew
"""

from matplotlib import use
use('agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
import datetime as dt
from os import system as sys
from temperature import TemperatureMap

years = range(2010, 2015)
months = range(1, 13)
days = [1, 15]
failed = []
dates = []
means = []
mins = []
maxes = []
stds = []

for year in years:
    for month in months:
        for day in days:
            try:
                date = dt.datetime(year, month, day)
                thismap = TemperatureMap(date)
                thismap.save()
                thismap.convert_scale('linear')
                fig = plt.figure(figsize=(16,12))
                #thismap.plot(cmap='gist_heat', vmin=0.6, vmax=2.5)#vmin=5.8, vmax=6.6)
                thismap.plot(cmap='RdYlBu_r', vmin=0.6, vmax=2.5)
                plt.title('Coronal temperature in MK')
                plt.colorbar()
                mapdir = '/media/huw/temperature_maps/1pars/maps/'
                error = sys('touch {}{:%Y/%m/%d/}'.format(mapdir, date))
                if error != 0:
                    sys('{0}{1:%Y}; {0}{1:%Y/%m}; {0}{1:%Y/%m/%d}'.format(
                            'mkdir '+mapdir, date))
                plt.savefig(mapdir+'{:%Y/%m/%d/%Y-%m-%dT%H:%M:%S}'.format(date))
                plt.close()
                dates.append(str(date))
                means.append(np.nanmean(thismap.data, dtype='float64'))
                mins.append(thismap.min())
                maxes.append(thismap.max())
                stds.append(np.nanstd(thismap.data, dtype='float64'))
            except:
                print 'Failed', date
                failed.append(date)

print 'Failed to create temperature maps for the following dates:'
for f in failed:
    print '\t', str(f)

fig = plt.figure(figsize=(18, 14))
ax = fig.add_subplot(1, 1, 1)
dates = mdates.datestr2num(dates)
plt.plot(dates, mins, color='blue', label='Minimum temperature')
plt.errorbar(dates, means, stds, color='black', label='Standard deviation')
plt.plot(dates, means, color='yellow', label='Mean temperature')
plt.plot(dates, maxes, color='red', label='Maximum temperature')
ax.xaxis_date()
fig.autofmt_xdate()
plt.legend(loc=4, fontsize=16)
plt.xlabel('Date', fontsize=24)
plt.ylabel('log(T)', fontsize=24)
plt.savefig('long_term_temp_graph')
plt.close()