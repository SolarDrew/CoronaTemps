# -*- coding: utf-8 -*-
"""
Created on Thu Jun  5 15:15:09 2014

@author: drew
"""

from __future__ import division
from matplotlib import use
use('agg')
import matplotlib.pyplot as plt
from matplotlib import cm, _cm
from matplotlib import patches
import numpy as np
import sunpy
from sunpy.map import Map, GenericMap
from sunpy.instr.aia import aiaprep
from scipy.io.idl import readsav as read
from sys import argv
from os import path, system, makedirs
import datetime as dt
from sunpy.time.timerange import TimeRange as tr
import glob
from itertools import product
from mpi4py import MPI
try:
    from fits import calc_fits
    print 'Fortran extension imported successfully'
except ImportError:
    print 'Current extension is broken, missing or incompatible.\n'\
        +'Compiling Fortran extension.'
    system(path.expanduser('f2py -c -m fits ~/CoronaTemps/fitsmodule.f90'))
    from fits import calc_fits


home = path.expanduser('~')
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()


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
    tresp = read('{}/CoronaTemps/aia_tresp'.format(home))
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

args = []
for a in argv[1:]:
    for f in [eval, sunpy.time.parse_time]:
        try:
            a = f(a)
            break
        except:
            continue
    args.append(a)

date, n_params, data_dir, datfile, submap, verbose, force_temp_scan = args

wlens = ['094', '131', '171', '193', '211', '335']
t0 = 5.6
thiswlen = None

if rank == 0:
    if datfile:
        images = {}
        f = open(datfile)
        # Loop through wavelengths
        for line in f:
            if line[:3] in wlens:
                allwlenmaps = []
                thiswlen = line[:3]
                print 'Loading {} files'.format(thiswlen)
            elif 'fits' in line:
                thismap = aiaprep(Map(line[:-1]))
                thismap.data /= thismap.exposure_time
                allwlenmaps.append(thismap)
            elif line.strip() in ['', '\n']:
                if thiswlen:
                    wlenmap = allwlenmaps[-1]
                    for thismap in allwlenmaps[:-1]:
                        wlenmap.data += thismap.data
                    wlenmap.data /= len(allwlenmaps)
                    images[thiswlen] = wlenmap
 
        images = [images[w] for w in wlens]
    else:
        images = []
        imagefiles = []
        for wl, wlen in enumerate(wlens):
            fits_dir = path.join(data_dir, '{:%Y/*/*}/{}'.format(date, wlen))
            if verbose: print 'Searching {} for AIA data'.format(fits_dir)
            timerange = tr(date - dt.timedelta(seconds=5),
                           date + dt.timedelta(seconds=11))
            ntimes = int(timerange.seconds())
            times = [time.start() for time in timerange.split(ntimes)]
            for time in times:
                filename = path.join(fits_dir,
                    'aia*{0:%Y?%m?%d}?{0:%H?%M?%S}*lev1?fits'.format(time))
                filelist = glob.glob(filename)
                if filelist != []:
                    if verbose: print 'File found: ', filelist[0]
                    imagefiles.append(filelist[0])
                    temp_im = aiaprep(Map(filelist[0]))
                    if submap:
                        temp_im = temp_im.submap(*submap)
                    temp_im.data /= temp_im.exposure_time # Can probably increase speed a bit by making this * (1.0/exp_time)
                    images.append(temp_im)
                    break
                else:
                    pass
            if len(images) < wl+1:
                if verbose: print 'Insufficient raw data - only found {} of 6 files'.format(len(images))

    # Normalise images to 171A if only using one parameter
    if n_params == 1:
        normim = images[2].data.copy()
        if verbose: print 'Normalising images'
        for i in range(len(wlens)):
            images[i].data /= normim
    header = images[2].meta.copy()
    images = np.array([im.data for im in images])

# Scatter image data to each process
if rank == 0:
    #[images[..., (p/size)*images.shape[2]:((p+1)/size)*images.shape[2]] \
    #    for p in range(size)]
    temp = []
    for p in range(size):
        mini = (p/size)*images.shape[2]
        maxi = ((p+1)/size)*images.shape[2]
        temp.append(images[..., mini:maxi])
        if verbose: print p, mini, maxi, images[..., mini:maxi].shape
    images = temp
    if verbose: print len(images), images[0].shape
else:
    images = None
images = comm.scatter(images, root=0)

# Get dimensions of image
x, y = images[0].shape
if verbose:
    print 'Image size, rank {}:'.format(rank), x, y
    print 'Image maxes, rank {}:'.format(rank), [im.max() for im in images]
n_wlens = images.shape[0]
temp = np.arange(t0, 7.01, 0.01)
if n_params == 1:
    # Assume a width of the gaussian DEM distribution and normalise the height
    widths = [0.1]
    heights = [1.0]
    parvals = temp
else:
    widths = np.arange(0.1, 0.8, 0.1)
    heights = [10.0**h for h in np.arange(20, 35, 2)]
    # TODO: check if either of the above are sensible ranges of numbers
    parvals = np.array([i for i in product(temp, widths, heights)])
n_vals = len(temp) * len(widths) * len(heights)

if rank == 0:
    try:
        if force_temp_scan:
            raise IOError
        model = np.memmap(filename='synth_emiss_{}pars'.format(n_params),
                          dtype='float32', mode='r', shape=(n_vals, n_wlens))
    except IOError:
        if verbose: print 'No synthetic emission data found. Re-scanning temperature range.'
        resp = load_temp_responses()
        if verbose: print resp.max(axis=1)
        logt = np.arange(0, 15.05, 0.05)
        delta_t = logt[1] - logt[0]
        model = np.memmap(filename='synth_emiss_{}pars'.format(n_params),
                          dtype='float32', mode='w+', shape=(n_vals, n_wlens))
        for p, params in enumerate(parvals):
            dem = gaussian(logt, *params)
            f = resp * dem
            model[p, :] = np.sum(f, axis=1) * delta_t
        if n_params == 1:
            normmod = model[:, 2]
            model /= normmod
        model.flush()
        if verbose: print model.max(axis=0)
else:
    model = None

model = comm.bcast(model, root=0)

if verbose:
    if rank == 0: print 'Calculating temperature values...',
    print rank, images.shape, model.shape, parvals.shape, n_vals, n_wlens, x, y, n_params
temps, fits = calc_fits(images, model, parvals, n_vals, n_wlens, x, y, n_params)
if verbose: print 'Done.'

# Get data all back in one place and save it
temps = comm.gather(temps, root=0)
if rank == 0:
    if verbose: print len(temps), temps[0].shape
    temp = np.zeros(shape=(x, y*size, n_params))
    for p in range(size):
        mini = (p/size)*temp.shape[1]
        maxi = ((p+1)/size)*temp.shape[1]
        temp[:, mini:maxi, :] = temps[p]
        if verbose: print p, mini, maxi, images[:, mini:maxi, :].shape
    temps = temp
    if verbose: print temps.shape
    tempmap = GenericMap(temps, header)
    tempmap.save(path.expanduser('~/CoronaTemps/temporary.fits'))
