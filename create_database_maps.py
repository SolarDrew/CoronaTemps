# -*- coding: utf-8 -*-
"""
Created on Wed Sep  3 16:13:05 2014

@author: drew
"""

from matplotlib import use
use('agg')
import matplotlib.pyplot as plt
from temperature import TemperatureMap
import datetime as dt
from os import system as sys

years = range(2010, 2015)
months = range(1, 13)
failed = []

for year in years:
    for month in months:
        try:
            date = dt.datetime(year, month, 1)
            thismap = TemperatureMap(date)
            thismap.save()
            thismap.convert_scale('linear')
            fig = plt.figure(figsize=(16,12))
            thismap.plot(cmap='gist_heat', vmin=0.6, vmax=2.5)#vmin=5.8, vmax=6.6)
            plt.title('Coronal temperature in MK')
            plt.colorbar()
            mapdir = '/media/huw/temperature_maps/1pars/maps/'
            error = sys('touch {}{:%Y/%m/%d/}'.format(mapdir, date))
            if error != 0:
                sys('{0}{1:%Y}; {0}{1:%Y/%m}; {0}{1:%Y/%m/%d}'.format(
                        'mkdir '+mapdir, date))
            plt.savefig(mapdir+'{:%Y/%m/%d/%Y-%m-%dT%H:%M:%S}'.format(date))
            plt.close()
        except:
            print 'Failed', date
            failed.append(date)

print 'Failed to create temperature maps for the following dates:'
for f in failed:
    print '\t', str(f)