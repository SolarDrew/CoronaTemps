# -*- coding: utf-8 -*-
"""
Created on Thu Jun  5 15:15:09 2014

@author: drew
"""

from matplotlib import use
use('agg')
import matplotlib.pyplot as plt
from matplotlib import cm, _cm
from matplotlib import patches
import numpy as np
import sunpy
from sunpy.map import Map, GenericMap
from sunpy.instr.aia import aiaprep
from sys import argv
from os import path, makedirs
import subprocess32 as subp
from scipy.io.idl import readsav as read
import glob


home = path.expanduser('~/macrospicules_paper/')
cortemps = path.join(home, 'CoronaTemps')


class TemperatureMap(GenericMap):
    def __init__(self, date=None, n_params=1, data_dir=None, maps_dir=None, 
                 fname=None, infofile=None, submap=None, verbose=False,
                 force_temp_scan=False):
        if (not fname and not date) or (fname and date):
            print """You must specify either a date and time for which to create
                temperatures or the name of a file containing a valid 
                TemperatureMap object."""
            return

        if date:
            date = sunpy.time.parse_time(date)
        
            if data_dir is None:
                data_dir = '/media/huw/SDO_data/'
            if maps_dir is None:
                maps_dir='/media/huw/temperature_maps/{}pars/'.format(n_params)
            
            fname = path.join(maps_dir, '{:%Y-%m-%dT%H_%M_%S}.fits'.format(date))

        if infofile:
            data_dir = None
            maps_dir = open(infofile).readline()[:-1]
            fname = path.join(maps_dir, '{:%Y-%m-%dT%H:%M:%S}.fits'.format(date))
            fname.replace('/images/', '/data/')

        if n_params != 1:
            fname = fname.replace('.fits', '_full.fits')

        if fname and not date:
            data_dir = path.dirname(fname)

        if verbose: print fname, data_dir

        try:
            newmap = Map(fname)
            GenericMap.__init__(self, newmap.data[..., 0], newmap.meta)
            self.goodness_of_fit = newmap.data[..., -1]
            if newmap.data.shape[2] != 2:
                self.dem_width = newmap.data[..., 1]
                self.emission_measure = newmap.data[..., 2]
        except ValueError:
            cmdargs = ["python", path.join(cortemps, 'create_tempmap.py'),
                date, n_params, data_dir, infofile, submap, verbose, force_temp_scan]
            if n_params != 1:
                cmdargs = ["mpiexec", "-n", 16] + cmdargs
            cmdargs = [str(cmd) for cmd in cmdargs]
            status = subp.call(cmdargs)
            newmap = Map(path.join(cortemps, 'temporary.fits'))
            subp.call(["rm", path.join(cortemps, 'temporary.fits')])
            data, meta = newmap.data, newmap.meta
            if verbose: print data.shape
            GenericMap.__init__(self, data[..., 0], meta)
            if data.shape[2] != 2:
                data[data == 0] = np.nan
                self.dem_width = data[..., 1]
                self.emission_measure = data[..., 2]
            self.goodness_of_fit = data[..., -1]
            if verbose:
                print self.shape
                print self.goodness_of_fit.shape
                if n_params != 1:
                    print self.dem_width.shape
                    print self.emission_measure.shape
            lowx, highx = (self.xrange[0] / self.scale['x'],
                           self.xrange[1] / self.scale['x'])
            lowy, highy = (self.yrange[0] / self.scale['y'],
                           self.yrange[1] / self.scale['y'])
            x_grid, y_grid = np.mgrid[lowx:highx-1, lowy:highy-1]
            r_grid = np.sqrt((x_grid ** 2.0) + (y_grid ** 2.0))
            outer_rad = (self.rsun_arcseconds * 1.5) / self.scale['x']
            self.data[r_grid > outer_rad] = None
            self.meta['date-obs'] = str(date)

        tmapcubehelix = _cm.cubehelix(s=2.8, r=0.7, h=2.0, gamma=1.0)
        cm.register_cmap(name='temphelix', data=tmapcubehelix)
        self.cmap = cm.get_cmap('temphelix')

        self.data_dir = data_dir
        self.maps_dir = maps_dir
        self.temperature_scale = 'log'
        self.region = None
        self.region_coordinate = {'x': 0.0, 'y': 0.0}
        if n_params == 3:
            self.n_params = 3
        else:
            self.n_params = 1

        return
    
    @classmethod
    def is_datasource_for(cls, data, header, **kwargs):
        return header.get('instrume', '').startswith('temperature')
    
    def region_map(self, region, mapsize=300, *args, **kwargs):
        """
        A function to take as input a hek record or similar and create a submap
        showing just the corresponding region
        """
        x, y = region['hpc_coord']
        newmap = self.submap([x-mapsize, x+mapsize], [y-mapsize, y+mapsize],
                            *args, **kwargs)
        
        self.region_coordinate = {'x': x, 'y': y}
        self.region = region

        return newmap
    
    def select_temps(self, mintemp, maxtemp):
        """
        Function to highlight user-defined temperatures
        """
        newdata = np.ones(self.data.shape) * np.NaN
        indices = np.where((self.data > mintemp) * (self.data < maxtemp))
        newdata[indices] = self.data[indices]
        
        return Map(newdata, self.meta.copy())
    
    def convert_scale(self, scale='linear'):
        if self.temperature_scale == scale:
            print "Temperatures are already measured on a {} scale.".format(
                scale)
            return
        elif scale == 'linear':
            self.data = (10.0 ** self.data) / 1.0e6
        elif scale == 'log':
            self.data = np.log10(self.data)
        
        self.temperature_scale = scale
        return
    
    def compare(self, display_wlen='171', context_wlen=None, extra_maps=[]):
        valid_wlens = ['94', '131', '171', '195', '211', '335', '304', 'hmi']
        if display_wlen.lower() not in valid_wlens:
            print "Display wavelength provided invalid or None."
            output = self.plot()#*temp_args, **temp_kwargs)
            return output
        save_output = True
        data_dir = self.data_dir
        maps_dir = self.maps_dir
        
        date = self.date
        nmaps = 2 + len(extra_maps)
        if context_wlen:
            nrows = 2
        else:
            nrows = 1
        
        fig = plt.figure(figsize=(24, 14))
        
        fig.add_subplot(nrows, nmaps, nmaps, axisbg='k')
        self.plot()#*temp_args, **temp_kwargs)
        plt.colorbar(orientation='horizontal')
        
        displaymap = Map(data_dir+'{0}/{1:%Y/%m/%d}/aia*{0}*t{1:%H?%M}*lev1?fits'\
            .format(display_wlen, date))
        if isinstance(displaymap, list):
            displaymap = displaymap[0]
        displaymap = aiaprep(displaymap)
        displaymap /= displaymap.exposure_time
        
        fig.add_subplot(nrows, nmaps, 1, axisbg='k')
        displaymap.plot()#*wlen_args, **wlen_kwargs)
        plt.colorbar(orientation='horizontal')
        
        if context_wlen and self.region != None:
            context_plot = fig.add_subplot(nrows, 1, nrows)
            contextmap = Map(data_dir+'{0}/{1:%Y/%m/%d}/aia*{0}*t{1:%H?%M}*lev1?fits'.format(context_wlen, date))
            if isinstance(contextmap, list):
                contextmap = contextmap[0]
            x, y = self.region_coordinate['x'], self.region_coordinate['y']
            contextmap = contextmap.submap([-1000, 1000], [y-300, y+300])
            # Need to figure out how to get 'subimsize' from self. Use the default 150'' for now
            #rect = patches.Rectangle([x-subdx, y-subdx], subimsize[0], subimsize[1], color='white', fill=False)
            rect = patches.Rectangle([x-150, y-150], 300, 300, color='white',
                                     fill=False)
            contextmap.plot()#*ctxt_args, **ctxt_kwargs)
            context_plot.add_artist(rect)
        
        for m, thismap in extra_maps:
            fig.add_subplot(nrows, nmaps, 3+m)
            thismap.plot()#*extr_args, **extr_kwargs)
        
        if save_output:
            savedir = path.join(maps_dir, 'maps/{:%Y/%m/%d}'.format(date))
            if not path.exists(savedir):
                makedirs(savedir)
            filename = path.join(maps_dir, '{%Y-%m-%dT%H:%M:%S}_with{}'.format(date, display_wlen))
            plt.savefig(filename)
            if self.region != None:
                reg_dir = path.join(maps_dir,
                                    'maps/region_maps/{}/'. format(self.region))
                if not path.exists(reg_dir):
                    makedirs(reg_dir)
                plt.savefig(path.join(reg_dir, '{:%Y-%m-%dT%H:%M:%S}'.format(date)))
            plt.close()
        else:
            plt.show()

        return
    
    def plot(self, vmin=None, vmax=None, *args, **kwargs):
        mean = np.nanmean(self.data, dtype=np.float64)
        std = np.nanstd(self.data, dtype=np.float64)
        if vmin is None:
            vmin = mean - (2.0 * std)
        if vmax is None:
            vmax = mean + (2.0 * std)
        
        GenericMap.plot(self, vmin=vmin, vmax=vmax, *args, **kwargs)
        
        return
    
    def save(self):
        date = sunpy.time.parse_time(self.date)
        if not path.exists(self.maps_dir):
            makedirs(self.maps_dir)
        fname = path.join(self.maps_dir,
                          '{:%Y-%m-%dT%H_%M_%S}.fits'.format(date))
        alldata = np.zeros((self.shape[0], self.shape[1], self.n_params+1))
        alldata[..., 0] = self.data
        if self.n_params != 1:
            fname = fname.replace('.fits', '_full.fits')
            alldata[..., 1] = self.dem_width
            alldata[..., 2] = self.emission_measure
        alldata[..., -1] = self.goodness_of_fit
        outmap = Map(alldata, self.meta.copy())
        outmap.save(fname, clobber=True)


    def min(self):
        return np.nanmin(self.data)

    def mean(self):
        return np.nanmean(self.data, dtype='float64')

    def max(self):
        return np.nanmax(self.data)

    def std(self):
        return np.nanstd(self.data, dtype='float64')

    def calculate_em(self, wlen='171', dz=100, model=False):
        """
        Calculate an approximation of the coronal EmissionMeasure using a given
        TemperatureMap object and a particular AIA channel.
    
        Parameters
        ----------
        tmap : CoronaTemps.temperature.TemperatureMap
            A TemperatureMap instance containing coronal temperature data
        wlen : {'94' | '131' | '171' | '193' | '211' | '335'}
            AIA wavelength used to approximate the emission measure. '171', '193'
            and '211' are most likely to provide reliable results. Use of other
            channels is not recommended.
        """
        # Load the appropriate temperature response function
        tresp = read('/imaps/holly/home/ajl7/CoronaTemps/aia_tresp')
        resp = tresp['resp{}'.format(wlen)]
    
        # Get some information from the TemperatureMap and set up filenames, etc
        tempdata = self.data.copy()
        tempdata[np.isnan(tempdata)] = 0.0
        date = sunpy.time.parse_time(self.date)
        if not model:
            data_dir = self.data_dir
            fits_dir = path.join(data_dir, '{:%Y/%m/%d}/{}'.format(date, wlen))
            filename = path.join(fits_dir,
                                 '*{0:%Y?%m?%d}?{0:%H?%M}*fits'.format(date))
            if wlen == '94': filename = filename.replace('94', '094')
    
            # Load and appropriately process AIA data
            filelist = glob.glob(filename)
            if filelist == []:
                print 'AIA data not found :('
                return
            aiamap = Map(filename)
            aiamap.data /= aiamap.exposure_time
            aiamap = aiaprep(aiamap)
            aiamap = aiamap.submap(self.xrange, self.yrange)
        else:
            fname = '/imaps/holly/home/ajl7/CoronaTemps/data/synthetic/{}/model.fits'.format(wlen)
            if wlen == '94': fname = fname.replace('94', '094')
            aiamap = Map(fname)

        # Create new Map and put EM values in it
        emmap = Map(self.data.copy(), self.meta.copy())
        indices = np.round((tempdata - 4.0) / 0.05).astype(int)
        indices[indices < 0] = 0
        indices[indices > 100] = 100
        #print emmap.shape, indices.shape, tempdata.shape, aiamap.shape, resp.shape
        emmap.data = np.log10(aiamap.data / resp[indices])
        #emmap.data = aiamap.data / resp[indices]

        emmapcubehelix = _cm.cubehelix(s=2.8, r=-0.7, h=1.4, gamma=1.0)
        cm.register_cmap(name='emhelix', data=emmapcubehelix)
        emmap.cmap = cm.get_cmap('emhelix')
    
        return emmap


sunpy.map.Map.register(TemperatureMap, TemperatureMap.is_datasource_for)

if __name__ == "__main__":
    date = sunpy.time.parse_time(argv[1])
    infofile = argv[2]    
    tmap = TemperatureMap(date, infofile=infofile)
    tmap.save()
    
    image_dir = open(infofile).readline()[:-1]
    fname = path.join(image_dir, '{:%Y-%m-%dT%H_%M_%S}'.format(date))
    print "Temperature map image saved to: {}".format(fname)
    
    fig = plt.figure(16, 12)
    tmap.plot()
    plt.colorbar(orientation='vertical')
    plt.savefig(fname)
    plt.close()
