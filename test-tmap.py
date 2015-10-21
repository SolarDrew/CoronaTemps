from matplotlib import use, rc
use('agg')
import matplotlib.pyplot as plt
from temperature import TemperatureMap as tm
from sunpy.map import Map
from sunpy.time import parse_time as pt
import numpy as np
from os.path import expanduser

date = '2014-01-11 13:00'

tmap = tm(date, n_params=1, verbose=True,
          data_dir=expanduser('~/SDO-data/'),
          maps_dir=expanduser('~/CoronaTemps/'))
#tmap.save()
