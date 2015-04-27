from matplotlib import use
use('agg')
import matplotlib.pyplot as plt
from os import path
from temperature import TemperatureMap as tmap


data_dir = path.expanduser('~/CoronaTemps/')
maps_dir = path.expanduser('~/CoronaTemps/')
thismap = tmap('2011-02-15', data_dir=data_dir, maps_dir=maps_dir, n_params=3, verbose=True)

print thismap.shape
print thismap.min(), thismap.mean(), thismap.max()

fig = plt.figure(figsize=(14, 14))
thismap.plot()
plt.savefig(path.expanduser('~/CoronaTemps/fullgausstemps'))
plt.close()
