from matplotlib import use
use('agg')
import matplotlib.pyplot as plt
import numpy as np
from os import path
from temperature import TemperatureMap as tmap
from sunpy import map


data_dir = path.expanduser('~/CoronaTemps/')
maps_dir = path.expanduser('~/CoronaTemps/')
data_dir = '/imaps/sspfs/archive/sdo/aia/fulldisk/data/'
thismap1 = tmap('2011-02-15', data_dir=data_dir, maps_dir=maps_dir, verbose=True)
thismap1.save()
thismap3 = tmap('2011-02-15', data_dir=data_dir, maps_dir=maps_dir, n_params=3,
                verbose=True, force_temp_scan=True)
#thismap3.save()
widths = map.Map(thismap3.dem_width, thismap3.meta)
EMs = map.Map(thismap3.emission_measure, thismap3.meta)

print thismap3.shape
print thismap3.dem_width.shape, np.nanmin(thismap3.dem_width), np.nanmean(thismap3.dem_width), np.nanmax(thismap3.dem_width)
print thismap3.emission_measure.shape, np.nanmin(thismap3.emission_measure), np.nanmean(thismap3.emission_measure), np.nanmax(thismap3.emission_measure)
print np.nanmin(thismap1.data), np.nanmean(thismap1.data), np.nanmax(thismap1.data)
print np.nanmin(thismap3.data), np.nanmean(thismap3.data), np.nanmax(thismap3.data)

fig = plt.figure(figsize=(14, 14))
thismap1.plot()
plt.colorbar()
plt.savefig(path.expanduser('~/CoronaTemps/1paramtemps'))
plt.close()

fig, ax = plt.subplots(1, 3, figsize=(42, 14))
plt.sca(ax[0])
thismap3.plot()
plt.colorbar()
plt.sca(ax[1])
widths.plot()#vmin=0.1, vmax=0.7)
plt.colorbar()
plt.sca(ax[2])
EMs.plot()#vmin=19, vmax=27)
plt.colorbar()
plt.savefig(path.expanduser('~/CoronaTemps/fullgausstemps'))
plt.close()
