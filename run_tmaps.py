from matplotlib import use
use('agg')
import sys
from os.path import join, expanduser
CThome = expanduser(join('~', 'macrospicules_paper', 'CoronaTemps'))
sys.path.append(CThome)
from temperature import TemperatureMap as TMap
from matplotlib import pyplot as plt
from sunpy.time import parse_time as parse

dates = ['2012-06-21 07:20:09']

for date in dates:
    data_path = join('/fastdata', 'sm1ajl', 'macrospicules',
        'aia*{w}a_{d.year:04d}_{d.month:02d}_{d.day:02d}?{d.hour:02d}_{d.minute:02d}_{d.second:02d}*.fits')
    map_path = join(CThome)# '{:%Y-%m-%dT%H_%M_%S}.fits'.format(date))
    print date
    print data_path
    print map_path
    thismap = TMap(date, data_path=data_path, map_path=map_path, verbose=True)
    thismap.save()
    #thismap = thismap.submap([-1200, -200], [-500, 500], units='data')
    fig = plt.figure(figsize=(32, 24))
    thismap.plot()
    plt.colorbar()
    plt.savefig('t_{:%Y-%m-%d_%H-%M-%S}'.format(parse(thismap.date)))
    plt.close()
