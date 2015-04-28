from matplotlib import use
use('agg')
import matplotlib.pyplot as plt
import numpy as np
from os import path
from temperature import TemperatureMap as tmap


data_dir = path.expanduser('~/CoronaTemps/')
maps_dir = path.expanduser('~/CoronaTemps/')
data_dir = '/imaps/sspfs/archive/sdo/aia/fulldisk/data/'
thismap1 = tmap('2011-02-15', data_dir=data_dir, maps_dir=maps_dir, verbose=True)
thismap1.save()
thismap3 = tmap('2011-02-15', data_dir=data_dir, maps_dir=maps_dir, n_params=3, verbose=True)
thismap3.save()

print thismap3.shape
print thismap3.gaussian_parameters
print np.nanmin(thismap1.data), np.nanmean(thismap1.data), np.nanmax(thismap1.data)
print np.nanmin(thismap3.data), np.nanmean(thismap3.data), np.nanmax(thismap3.data)

fig = plt.figure(figsize=(14, 14))
thismap1.plot()
plt.colorbar()
plt.savefig(path.expanduser('~/CoronaTemps/1paramtemps'))
plt.close()

fig = plt.figure(figsize=(14, 14))
thismap3.plot()
plt.colorbar()
plt.savefig(path.expanduser('~/CoronaTemps/fullgausstemps'))
plt.close()
