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

years = range(2010, 2015)
months = range(1, 13)

for year in years:
    for month in months:
        date = dt.datetime(year, month, 1)
        thismap = TemperatureMap(date)
        thismap.save()
        thismap.convert_scale('linear')
        fig = plt.figure(figsize=(16,12))
        thismap.plot(cmap='gist_heat', vmin=0.6, vmax=2.5)#vmin=5.8, vmax=6.6)
        plt.title('Coronal temperature in MK')
        plt.colorbar()
        plt.savefig(date)
        plt.close()