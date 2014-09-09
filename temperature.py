# -*- coding: utf-8 -*-
"""
Created on Thu Jun  5 15:15:09 2014

@author: drew
"""

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import patches
import numpy as np
import datetime as dt
import sunpy
from sunpy.net import vso
from sunpy.map import Map, GenericMap
from sunpy.instr.aia import aiaprep
from scipy.io.idl import readsav as read
from os import system as sys
from astropy import units as u
#import numexpr as ne
try:
    from fits import calc_fits
except ImportError:
    print 'Current extension is broken, missing or incompatible.\n'\
        +'Compiling Fortran extension.'
    sys('f2py -c -m fits fitsmodule.f90')
    from fits import calc_fits

home = '/media/huw/'

def gaussian(x, mean=0.0, std=1.0, amp=1.0):
    """Simple function to return a Gaussian distribution"""
    if isinstance(x, list):
        x = np.array(x)
    power = -((x - mean) ** 2.0) / (2.0 * (std ** 2.0))
    f = amp * np.exp(power)
    if amp == 1:
        f = f / max(f)
    return f


def load_temp_responses(n_wlens=6, corrections=True):
    resp = np.zeros((n_wlens, 301))
    try:
        tresp = read(home + 'aia_tresp')
    except IOError:
        tresp = read('/imaps/holly/home/ajl7/aia_tresp')
    resp[0, 80:181] = tresp['resp94']
    resp[1, 80:181] = tresp['resp131']
    resp[2, 80:181] = tresp['resp171']
    resp[3, 80:181] = tresp['resp193']
    resp[4, 80:181] = tresp['resp211']
    resp[5, 80:181] = tresp['resp335']
    if n_wlens > 6:
        resp[6, 80:181] = tresp['resp304']
    if corrections:
        # Add empirical correction factor for 9.4nm response function below log(T)=6.3
        # (see Aschwanden et al 2011)
        resp[0:126, 0] = resp[0:126, 0]*6.7
    
    return resp


def find_temp(images, t0=5.6, force_temp_scan=False, maps_dir=home+'temperature_maps/'):
    x, y = images[0].shape
    n_wlens = len(images)
    n_temps = int((7.0 - t0) / 0.01) + 1
    temp = np.arange(t0, 7.01, 0.01)
    
    try:
        if force_temp_scan:
            raise IOError
        model = np.memmap(filename=home+'synth_emiss_1pars', dtype='float32',
                          mode='r', shape=(n_temps, n_wlens))
    except IOError:
        print 'No synthetic emission data found. Re-scanning temperature range.'
        resp = load_temp_responses()
        logt = np.arange(0, 15.05, 0.05)
        # Assume a width of the gaussian DEM distribution and normalise the height
        width = 0.1
        height = 1.0
        delta_t = logt[1] - logt[0]
        model = np.memmap(filename=home+'synth_emiss_1pars', dtype='float32',
                          mode='w+', shape=(n_temps, n_wlens))
        for t, meantemp in enumerate(temp):
            dem = gaussian(logt, meantemp, width, height)
            f = resp * dem
            model[t, :] = np.sum(f, axis=1) * delta_t ### CHECK THIS AXIS!
            normmod = model[t, 2]
            model[t, :] = model[t, :] / normmod
        model.flush()
    ims_array = np.array([im.data for im in images])
    print 'Calculating temperature values...',
    temps, fits = calc_fits(ims_array, model, temp, n_temps, n_wlens, x, y, 1)
    print 'Done.'
    tempmap = temps[:, :, 0], images[2].meta.copy(), fits
    # TODO: figure out how to change things in the header and save them.
    
    return tempmap


"""def find_temp_3params(aiamaps, t0, force_temp_scan=False, maps_dir=home+'temperature_maps/'):
    """"""from mpi4py import MPI
    #fcomm = MPI.COMM_WORLD.py2f()
    
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()""""""
    
    if rank == 0:
        ims_array = np.array([im.data for im in aiamaps])
        n = ims_array.shape[-1]
        print ims_array.shape, n
        images = np.array([ims_array[:, :, r*(n/size):(r+1)*(n/size)] for r in range(size)])
    images = comm.scatter(images, root=0)
    print rank, images.shape
    
    #x, y = images[0].shape
    x, y = images.shape[-2], images.shape[-1]
    n_wlens = images.shape[0]#len(images)
    temp = np.arange(t0, 7.01, 0.01)
    n_temps = len(temp)
    wid = np.arange(0.1, 1.1, 0.1) # Just copying Aschwanden's range here
    n_widths = len(wid)
    hei = np.arange(15, 30)
    n_heights = len(hei)
    n_vals = n_temps * n_widths * n_heights
    params = np.zeros((n_vals, 3))
    
    try:
        if force_temp_scan:
            raise IOError
        model = np.memmap(filename=home+'synth_emiss_3pars', dtype='float32',
                          mode='r',shape=(n_vals, 1, 1, n_wlens))
    except IOError:
        if rank == 0:
            resp = load_temp_responses()
            logt = np.arange(0, 15.05, 0.05)
            delta_t = logt[1] - logt[0]
            model = np.memmap(filename=home+'synth_emiss_3pars', dtype='float32',
                              mode='w+', shape=(n_vals, 1, 1, n_wlens))
            #   For given gaussian width, height and centre:
            for w, width in enumerate(wid):
                for h, height in enumerate(hei):
                    for t, meantemp in enumerate(temp):
                        # Define linear index in model for this set of parameters
                        i = w + h + t
                        # Store these parameters for use later
                        params[i, :] = [meantemp, height, width]
                        # For each channel k:
                        # intensity(k, params) <- sum(temperature_response(k, params) * DEM(params)) * dT
                        dem = gaussian(logt, meantemp, width, height)
                        model[i, 0, 0, :] = np.sum(resp * dem, axis=1) * delta_t
                        #normmod = model[t, 0, 0, 2]
                        #model[i, 0, 0, :] = model[i, 0, 0, :] / normmod
            model.flush()
        model = comm.bcast(model, root=0)
    # ----- Load raw AIA data -----
    #   if lvl 1.5 data is not present:
    #       if lvl 1 data is not present:
    #           download lvl 1 data
    #       process lvl 1 data to lvl 1.5
    #       save lvl 1.5 data
    #   load lvl 1.5 data (into MapCube?)
    # For now, stick with what I was doing before
    """"""ims_array = np.memmap(filename=home+'images', dtype='float32', mode='w+',
                      shape=(1, x, y, n_wlens))""""""
    #ims_array = np.array([im.data for im in images])
    """"""for i, im in enumerate(images):
        ims_array[0, :, :, i] = im.data""""""
    
    # Create MapCube containing separate maps for temperature, emission measure
    # and DEM width
    # Possibly also add a map for density
    
    print 'Finding best Gaussian parameter values...',
    #best_params = calc_fits(fcomm, ims_array, model, params, n_vals, n_wlens, x, y, 3)
    best_params = calc_fits(images, model, params, n_vals, n_wlens, x, y, 3)
    #best_params = calc_fits_py(ims_array, model, params, n_vals, n_wlens, x, y)
    print 'Done.'
    print rank, best_params.shape
    #tempmap = best_params[:, :, 0], images[2].meta.copy()
    #print tempmap[0].shape
    #print tempmap[0].min(), tempmap[0].mean(), tempmap[0].max()
    tempdata = best_params[:, :, 0]
    tempmap = comm.gather(tempdata, root=0), aiamaps[2]
    # TODO: figure out how to change things in the header and save them.
    
    """"""tempmap = comm.Gather(tempmap, root=0)
    if rank == 0:
        print images.shape""""""
    
    return tempmap"""


def create_tempmap(date, n_params=1, data_dir=home+'SDO_data/',
                   maps_dir=home+'temperature_maps/'):
    wlens = ['94', '131', '171', '193', '211', '335']
    t0 = 5.6
    images = []
    #imdates = {}
    
    print 'Finding data for {}.'.format(date.date())
    # Loop through wavelengths
    for wl, wlen in enumerate(wlens):
        #print 'Finding {}A data...'.format(wlen),
        fits_dir = data_dir + '{}/{:%Y/%m/%d}/'.format(wlen, date)
        filename = fits_dir + 'aia*{0}*{1:%Y?%m?%d}?{1:%H?%M}*lev1?fits'.format(wlen, date)
        temp_im = Map(filename)
        # Download data if not enough found
        client = vso.VSOClient()
        if temp_im == []:
            print 'File not found. Downloading from VSO...'
            # Wavelenth needs to be a quantity for VSO query
            wQuant = u.Quantity(value=int(wlen), unit='Angstrom')
            qr = client.query(vso.attrs.Time(date,# - dt.timedelta(seconds=6),
                                             date + dt.timedelta(seconds=12)),#6)),
                              vso.attrs.Wave(wQuant, wQuant),
                              vso.attrs.Instrument('aia'),
                              vso.attrs.Provider('JSOC'))
            res = client.get(qr, path=fits_dir+'{file}', site='NSO').wait()
            temp_im = Map(res)
        if temp_im == []:
            print 'Downloading failed.'
            print res, len(qr), qr
            return np.zeros((512, 512))
        if isinstance(temp_im, list):
            temp_im = temp_im[0]
        # TODO: save out level 1.5 data so it can be loaded quickly.
        temp_im = aiaprep(temp_im)
        temp_im.data = temp_im.data / temp_im.exposure_time # Can probably increase speed a bit by making this * (1.0/exp_time)
        images.append(temp_im)
        #imdates[wlen] = temp_im.date
    
    normim = images[2].data.copy()
    # Normalise images to 171A
    print 'Normalising images'
    for i in range(len(wlens)):
        images[i].data = images[i].data / normim
    
    # Produce temperature map
    if n_params == 1:
        tempmap = find_temp(images, t0)#, force_temp_scan=True)
    else:
        #tempmap = find_temp_3params(images, t0)
        pass

    return tempmap


def calculate_temperatures(date, n_params=1, data_dir=home+'SDO_data/',
                            maps_dir=home+'temperature_maps/', n_procs=4):
    wlens = ['94', '131', '171', '193', '211', '335']
    client = vso.VSOClient()
    print 'Finding data for {}.'.format(date.date())
    # Loop through wavelengths
    for wl, wlen in enumerate(wlens):
        #print 'Finding {}A data...'.format(wlen),
        fits_dir = data_dir + '{}/{:%Y/%m/%d}/'.format(wlen, date)
        filename = fits_dir + 'aia*{0}*{1:%Y?%m?%d}?{1:%H?%M}*lev1?fits'.format(wlen, date)
        temp_im = Map(filename)
        # Download data if not enough found
        if temp_im == []:
            print 'File not found. Downloading from VSO...'
            qr = client.query(vso.attrs.Time(date,# - dt.timedelta(seconds=6),
                                             date + dt.timedelta(seconds=12)),#6)),
                              vso.attrs.Wave(wlen, wlen),
                              vso.attrs.Instrument('aia'),
                              vso.attrs.Provider('JSOC'))
            res = client.get(qr, path=fits_dir+'{file}', site='NSO',
                             methods=['URL_FILE_Rice']).wait()

    n_wlens = len(wlens)
    temp = np.arange(5.6, 7.01, 0.01)
    n_temps = len(temp)
    if n_params == 1:
        wid = [0.1]
        hei = [1.0]
    else:
        wid = np.arange(0.1, 1.1, 0.1) # Just copying Aschwanden's range here
        hei = np.arange(15, 31)
    n_widths = len(wid)
    n_heights = len(hei)
    n_vals = n_temps * n_widths * n_heights
    
    try:
        #raise IOError
        model = np.memmap(filename=home+'synth_emiss_{}pars'.format(n_params),
                          dtype='float32', mode='r', shape=(n_vals, n_wlens))
        params = np.loadtxt(home+'gaussian_parameters_{}'.format(n_params))
    except IOError:
        resp = load_temp_responses()
        logt = np.arange(0, 15.05, 0.05)
        delta_t = logt[1] - logt[0]
        model = np.memmap(filename=home+'synth_emiss_{}pars'.format(n_params),
                          dtype='float32', mode='w+', shape=(n_vals, n_wlens))
        params = np.zeros((n_vals, n_params))
        # Define linear index in model for this set of parameters
        i = 0
        #   For given gaussian width, height and centre:
        for width in wid:
            for height in hei:
                for meantemp in temp:
                    # For each channel k:
                    # intensity(k, params) <- sum(temperature_response(k, params) * DEM(params)) * dT
                    dem = gaussian(logt, float(meantemp), width, height)
                    model[i, :] = np.sum(resp * dem, axis=1) * delta_t
                    normmod = model[i, 2]
                    model[i, :] = model[i, :] / normmod
                    # Store these parameters for use later
                    params[i, :] = [meantemp, width, height]
                    i += 1
        model.flush()
        np.savetxt(home+'gaussian_parameters_{}'.format(n_params), params)

    sys('mpirun -np {} -host localhost python parallel_ext.py {}T{} {} {}'.format(
            n_procs, date.date(), date.time(), n_params, n_vals))
    
    return


class TemperatureMap(GenericMap):
    def __init__(self, date=None, n_params=1, data_dir=None, maps_dir=None, 
                 fname=None):
        if (not fname and not date) or (fname and date):
            print """"You must specify either a date and time for which to create
                temperatures or the name of a file containing a valid 
                TemperatureMap object."""
            return

        if date:
            date = sunpy.time.parse_time(date)
        
            if data_dir is None:
                data_dir = '/media/huw/SDO_data/'
            if maps_dir is None:
                maps_dir='/media/huw/temperature_maps/{}pars/'.format(n_params)
            
            fname = maps_dir+'data/{:%Y/%m/%d/%Y-%m-%dT%H:%M:%S}.fits'.format(date)

        try:
            #raise ValueError
            #fname = maps_dir+'data/{:%Y/%m/%d/%Y-%m-%dT%H:%M:%S}.fits'.format(date)
            newmap = Map(fname)
            GenericMap.__init__(self, newmap.data, newmap.meta)
            #self.data = newmap.data
            #self.meta = newmap.meta
        except ValueError:
            if n_params == 3:
                calculate_temperatures(date, n_params, data_dir, maps_dir, 8)
                newmap = Map('temp.fits')
                GenericMap.__init__(self, newmap.data[:, :, 0], newmap.meta)
            else:
                #data, meta = create_tempmap(date, n_params, data_dir, maps_dir)
                data, meta, fit = create_tempmap(date, n_params, data_dir, maps_dir)
                GenericMap.__init__(self, data, meta)
                centre_x = self.reference_pixel['x']
                centre_y = self.reference_pixel['y']
                x_grid, y_grid = np.mgrid[-centre_x:centre_x-1, -centre_y:centre_y-1]
                r_grid = np.sqrt((x_grid ** 2.0) + (y_grid ** 2.0))
                self.data[r_grid > centre_x * 1.15] = np.NaN
                self.data[self.data == 0.0] = np.NaN
            #self.data = data
            #self.meta = meta

        #self.date = date
        self.meta['date-obs'] = str(date)
        self.data_dir = data_dir
        self.maps_dir = maps_dir
        self.temperature_scale = 'log'
        self.cmap = cm.coolwarm
        self.region = None
        self.region_coordinate = {'x': 0.0, 'y': 0.0}
        if n_params == 3:
            self.n_params = 3
        else:
            self.n_params = 1

        """if n_params != 1:
            alldata = self.data
            self.data = alldata[:, :, 0]
            self.EM = alldata[:, :, 1]
            self.width = alldata[:, :, 2]"""
        
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
        splitmap = TemperatureMap(np.ones(self.data.shape)*np.NaN,
                                  self.meta.copy())
        indices = np.where((self.data > mintemp) * (self.data < maxtemp))
        splitmap.data[indices] = splitmap.data[indices]
        
        return splitmap
    
    def convert_scale(self, scale='linear'):
        if self.temperature_scale == scale:
            print "Temperatures are already measured on a {} scale.".format(
                scale)
            return
        elif scale == 'linear':
            self.data = (10.0 ** self.data) / 1.0e6
        elif scale == 'log':
            self.data = np.log(self.data)
        
        self.temperature_scale = scale
        return
    
    def compare(self, display_wlen='171', context_wlen=None, extra_maps=[]):
        #        temp_args=None, temp_kwargs=None,
        #        wlen_args=None, wlen_kwargs=None,
        #        ctxt_args=None, ctxt_kwargs=None,
        #        extr_args=None, extr_kwargs=None):
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
            error = sys('touch '+maps_dir+'maps/{:%Y/%m/%d/} > shelloutput.txt'.format(date))
            if error != 0:
                sys('{0}{1:%Y}; {0}{1:%Y/%m}; {0}{1:%Y/%m/%d} > shelloutput.txt'\
                        .format('mkdir '+maps_dir+'maps/', date))
            filename = maps_dir+'maps/{:%Y/%m/%d/%Y-%m-%dT%H:%M:%S}_with{}'.\
                    format(date, display_wlen)
            plt.savefig(filename)
            if self.region != None:
                reg_dir = maps_dir + 'maps/region_maps'
                reg_dir = reg_dir + '/{}/'. format(self.region)
                error = sys('touch ' + reg_dir + ' > shelloutput.txt')
                if error != 0:
                    sys('mkdir ' + reg_dir + ' > shelloutput.txt')
                plt.savefig(reg_dir+'{:%Y-%m-%dT%H:%M:%S}'.format(date))
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
        error = sys('touch '+self.maps_dir+'data/{:%Y/%m/%d/} > shelloutput.txt'.format(date))
        if error != 0:
            sys('{0}{1:%Y}; {0}{1:%Y/%m}; {0}{1:%Y/%m/%d} > shelloutput.txt'.format(
                'mkdir '+self.maps_dir+'data/', date))
        GenericMap.save(self, self.maps_dir+'data/{:%Y/%m/%d/%Y-%m-%dT%H:%M:%S}.fits'.format(date), clobber=True)

sunpy.map.Map.register(TemperatureMap, TemperatureMap.is_datasource_for)
