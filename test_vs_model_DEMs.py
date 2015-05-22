# -*- coding: utf-8 -*-
"""
Script to produce synthetic AIA data based on arbitrary model DEMs and test the
results of the tempmap code against the model.

Created on Mon Jul 28 16:34:28 2014

@author: Drew Leonard
"""

import numpy as np
from matplotlib import use
use('agg')
import matplotlib.pyplot as plt
from matplotlib import patches
import sunpy
from sunpy.map import Map
from temperature import TemperatureMap
from utils import gaussian, load_temp_responses
from os import path, makedirs
import subprocess32 as subp
from itertools import product

# Decide whether to assess single-parameter or full-Gaussian method
n_pars = 3

# Define CoronaTemps home folder and output folder
CThome = path.join(path.expanduser('~'), 'CoronaTemps')
outdir = path.join(CThome, 'validation', '{}pars'.format(n_pars))
tmap_script = path.join(CThome, 'create_tempmap.py')
if not path.exists(outdir): makedirs(outdir)

# Define parameter ranges
temps = np.arange(4.6, 7.405, 0.005)
#temps = np.arange(5.6, 7.005, 0.005)
widths = np.array([0.01, 0.1, 0.5])#np.arange(0.01, 0.605, 0.005) # Just copying Aschwanden's range here
heights = 10 ** np.arange(18, 37, 0.05)
#heights = 10 ** np.arange(20, 35, 0.05)
print heights
n_temps = len(temps)
n_widths = len(widths)
n_heights = len(heights)
parvals = np.array([i for i in product(temps, widths, heights)])
print parvals.shape
n_vals = n_temps * n_widths * n_heights
print n_temps, n_widths, n_heights, n_vals, n_vals * 6

# Create model DEMs and synthetic emission
emission = np.zeros((6, n_temps, n_widths, n_heights))
print emission.shape
logt = np.arange(0, 15.05, 0.05)
resp = load_temp_responses()
delta_t = logt[1] - logt[0]
for p, params in enumerate(parvals):
    dem = gaussian(logt, *params)
    f = resp * dem
    t = np.where(temps == params[0])[0][0]
    w = np.where(widths == params[1])[0][0]
    h = np.where(heights == params[2])[0][0]
    emission[:, t, w, h] = np.sum(f, axis=1) * delta_t

#emission = emission / emission[2, :, :, :]

# Load AIA response functions
resp = load_temp_responses()

# Load unnessecary map for its metadata
voidmap = Map(sunpy.AIA_171_IMAGE)
mapmeta = voidmap.meta

#rect = patches.Rectangle([25.0, 5.6], 1.0, 1.0, color='black', fill=True, clip_on=False)
# Run synthetic data through 1param tempmap method
for w, wid in enumerate(widths):#heights):
    print 'Width:', wid
    fig = plt.figure(figsize=(30, 12))
    for wl, wlength  in enumerate(['94', '131', '171', '193', '211', '335']):
        #emiss = Map(emission[wl, :, :, w], mapmeta)
        emiss = Map(emission[wl, :, w, :], mapmeta)
        emiss.cmap = sunpy.cm.get_cmap('sdoaia{}'.format(wlength))
        emiss.meta['naxis1'] = emiss.shape[1]
        emiss.meta['naxis2'] = emiss.shape[0]
        #emiss.meta['cdelt1'] = widths[1] - widths[0]
        emiss.meta['cdelt1'] = np.log10(heights[1]) - np.log10(heights[0])
        emiss.meta['cdelt2'] = temps[1] - temps[0]
        #emiss.meta['crval1'] = widths[0]
        emiss.meta['crval1'] = np.log10(heights[0])
        emiss.meta['crval2'] = temps[0]
        emiss.meta['crpix1'] = 0.5
        emiss.meta['crpix2'] = 0.5
        if wlength == '94': wlength = '094'
        fits_dir = path.join(CThome, 'data', 'synthetic', wlength)
        if not path.exists(fits_dir): makedirs(fits_dir)
        emiss.save(path.join(fits_dir, 'model.fits'), clobber=True)
        ax = fig.add_subplot(1, 6, wl+1)
        emiss.plot(aspect='auto')#vmin=0.0, vmax=1.0)
        plt.title('{}'.format(wlength))
        plt.xlabel('Input EM')
        plt.ylabel('Input log(T)')
        plt.colorbar()
        #fig.gca().add_artist(rect)
        plt.axvline(20.0, color='white')
        plt.axvline(35.0, color='white')
        plt.axhline(5.6, color='white')
        plt.axhline(7.0, color='white')
    #plt.savefig(path.join(outdir, 'model_emission_h={}'.format(np.log10(wid)).replace('.', '_')))
    plt.savefig(path.join(outdir, 'model_emission_w={}'.format(wid).replace('.', '_')))
    plt.close()
    
    #images = [Map(emission[i, :, :, w], mapmeta) for i in range(6)]
    images = [Map(emission[i, :, w, :], mapmeta) for i in range(6)]
    cmdargs = "mpiexec -n 10 python {} model {} {} {} {} {} {}".format(
        tmap_script, n_pars, path.join(CThome, 'data'),
        None, None, True, True).split()
    status = subp.call(cmdargs)
    newmap = Map(path.join(CThome, 'temporary.fits'))
    subp.call(["rm", path.join(CThome, 'temporary.fits')])
    data, meta = newmap.data[..., 0], newmap.meta
    wdata = newmap.data[..., 1]
    emdata = newmap.data[..., 2]
    print np.max(newmap.data[..., -1]), np.max(np.log10(newmap.data[..., -1]))
    fitsmap = Map(np.log10(newmap.data[..., -1]), newmap.meta.copy())
    print fitsmap.max()
    newmap.data = data
    #print '\n-----------MINMAX: ', newmap.min(), newmap.max(), newmap.shape, 
    #print newmap.data[newmap.data == 0].shape, '----------\n'

    #truetemp = np.array(list(temps)*n_widths).reshape((n_widths, n_temps)).T
    truetemp = np.array(list(temps)*n_heights).reshape((n_heights, n_temps)).T
    truew = np.ones(shape=(n_temps, n_heights)) * wid
    #print 'truew', truew.min(), truew.mean(), truew.max()
    trueem = np.log10(np.array(list(heights)*n_temps).reshape(n_temps, n_heights))
    diff = Map((abs(truetemp - data) / truetemp) * 100, newmap.meta.copy())
    diffw = Map((abs(truew - wdata) / truew) * 100, newmap.meta.copy())
    #print 'wdata', wdata.min(), wdata.mean(), wdata.max()
    #print 'diffw', diffw.min(), diffw.mean(), diffw.max()
    diffem = Map((abs(trueem - emdata) / trueem) * 100, newmap.meta.copy())
    
    print wid, newmap.xrange, newmap.yrange, newmap.scale
    print wid, diff.xrange, diff.yrange, diff.scale
    fig = plt.figure(figsize=(24, 12))
    fig.add_subplot(1, 3, 1)
    newmap.plot(cmap='coolwarm', vmin=5.6, vmax=7.0, aspect='auto')
    plt.colorbar()
    plt.title('Solution log(T)', fontsize=28)
    plt.ylabel('Input log(T)', fontsize=24)
    plt.xlabel('Input EM', fontsize=24)#width', fontsize=24)
    #fig.gca().add_artist(rect)
    plt.axvline(20.0, color='white')
    plt.axvline(35.0, color='white')
    plt.axhline(5.6, color='white')
    plt.axhline(7.0, color='white')
    
    ax = fig.add_subplot(1, 3, 2)
    print 'diff', diff.min(), diff.max()
    print np.nanmin(diff.data), np.nanmax(diff.data)
    diff.plot(cmap='RdYlGn_r', aspect='auto')#, vmin=diff.min(), vmax=diff.max())
    plt.colorbar()
    plt.title('Difference from input (%)', fontsize=28)
    plt.xlabel('Input EM', fontsize=24)
    #fig.gca().add_artist(rect)
    plt.axvline(20.0, color='white')
    plt.axvline(35.0, color='white')
    plt.axhline(5.6, color='white')
    plt.axhline(7.0, color='white')
    
    ax = fig.add_subplot(1, 3, 3)
    print 'fits', fitsmap.min(), fitsmap.max()
    print np.nanmin(fitsmap.data), np.nanmax(fitsmap.data)
    fitsmap.plot(cmap='cubehelix', aspect='auto')
    plt.colorbar()
    plt.title('log(Goodness-of-fit)', fontsize=28)
    plt.xlabel('Input EM', fontsize=24)
    #fig.gca().add_artist(rect)
    plt.axvline(20.0, color='white')
    plt.axvline(35.0, color='white')
    plt.axhline(5.6, color='white')
    plt.axhline(7.0, color='white')
    #plt.savefig(path.join(outdir, 'tempsolutions_em={}'.format(np.log10(wid)).replace('.', '_')))
    plt.savefig(path.join(outdir, 'tempsolutions_wid={:.3f}'.format(wid).replace('.', '_')))
    plt.close()

    wdata = Map(wdata, newmap.meta.copy())
    print wid, wdata.xrange, wdata.yrange, wdata.scale
    print wid, diffw.xrange, diffw.yrange, diffw.scale
    fig = plt.figure(figsize=(24, 12))
    ax = fig.add_subplot(1, 3, 1)
    wdata.plot(cmap='coolwarm', vmin=widths[0], vmax=widths[-1], aspect='auto')
    #fig.gca().add_artist(rect)
    plt.axvline(20.0, color='white')
    plt.axvline(35.0, color='white')
    plt.axhline(5.6, color='white')
    plt.axhline(7.0, color='white')
    plt.colorbar()
    plt.title('Solution width', fontsize=28)
    plt.ylabel('Input log(T)', fontsize=24)
    plt.xlabel('Input EM', fontsize=24)#width', fontsize=24)
    
    ax = fig.add_subplot(1, 3, 2)
    print 'diffw', diffw.min(), diffw.max()
    print np.nanmin(diffw.data), np.nanmax(diffw.data)
    diffw.plot(cmap='RdYlGn_r', vmin=diffw.min(), vmax=diffw.max(), aspect='auto')
    #fig.gca().add_artist(rect)
    plt.axvline(20.0, color='white')
    plt.axvline(35.0, color='white')
    plt.axhline(5.6, color='white')
    plt.axhline(7.0, color='white')
    plt.colorbar()
    plt.title('Difference from input (%)', fontsize=28)
    plt.xlabel('Input EM', fontsize=24)

    ax = fig.add_subplot(1, 3, 3)
    fitsmap.plot(cmap='cubehelix', aspect='auto')
    #fig.gca().add_artist(rect)
    plt.axvline(20.0, color='white')
    plt.axvline(35.0, color='white')
    plt.axhline(5.6, color='white')
    plt.axhline(7.0, color='white')
    plt.colorbar()
    plt.title('log(Goodness-of-fit)', fontsize=28)
    plt.xlabel('Input EM', fontsize=24)
    
    plt.savefig(path.join(outdir, 'widsolutions_wid={:.3f}'.format(wid).replace('.', '_')))
    plt.close()

    emdata = Map(emdata, newmap.meta.copy())
    print wid, emdata.xrange, emdata.yrange, emdata.scale
    print wid, diffem.xrange, diffem.yrange, diffem.scale
    print wid, fitsmap.xrange, fitsmap.yrange, fitsmap.scale
    fig = plt.figure(figsize=(24, 12))
    ax = fig.add_subplot(1, 3, 1)
    print 'emdata', emdata.max(), emdata.min()
    print np.nanmin(emdata.data), np.nanmax(emdata.data)
    emdata.plot(cmap='coolwarm', aspect='auto',
                vmin=np.log10(heights[0]), vmax=np.log10(heights[-1]))
    #fig.gca().add_artist(rect)
    plt.axvline(20.0, color='white')
    plt.axvline(35.0, color='white')
    plt.axhline(5.6, color='white')
    plt.axhline(7.0, color='white')
    plt.colorbar()
    plt.title('Solution EM', fontsize=28)
    plt.ylabel('Input log(T)', fontsize=24)
    plt.xlabel('Input EM', fontsize=24)#width', fontsize=24)
    
    ax = fig.add_subplot(1, 3, 2)
    print 'diffem', diffem.min(), diffem.max()
    print np.nanmin(diffem.data), np.nanmax(diffem.data)
    diffem.plot(cmap='RdYlGn_r', aspect='auto',
                vmin=diffem.min(), vmax=diffem.max())
    #fig.gca().add_artist(rect)
    plt.axvline(20.0, color='white')
    plt.axvline(35.0, color='white')
    plt.axhline(5.6, color='white')
    plt.axhline(7.0, color='white')
    plt.colorbar()
    plt.title('Difference from input (%)', fontsize=28)
    plt.xlabel('Input EM', fontsize=24)
    
    ax = fig.add_subplot(1, 3, 3)
    fitsmap.plot(cmap='cubehelix', aspect='auto')
    #fig.gca().add_artist(rect)
    plt.axvline(20.0, color='white')
    plt.axvline(35.0, color='white')
    plt.axhline(5.6, color='white')
    plt.axhline(7.0, color='white')
    plt.colorbar()
    plt.title('log(Goodness-of-fit)', fontsize=28)
    plt.xlabel('Input EM', fontsize=24)

    plt.savefig(path.join(outdir, 'emsolutions_wid={:.3f}'.format(wid).replace('.', '_')))
    plt.close()

"""w = np.where((widths > 0.097)*(widths < 0.103))
dataslice = data[:, w].reshape(len(temps))
diffslice = diff[:, w].reshape(len(temps))
fitslice = fits[:, w].reshape(len(temps))

fig = plt.figure(figsize=(16, 12))
plt.plot(temps, dataslice)
plt.title('Solution log(T) at width=0.1', fontsize=28)
plt.xlabel('Input log(T)', fontsize=24)
plt.ylabel('Solution log(T)', fontsize=24)
plt.savefig('/home/drew/Dropbox/euroscipy/dataslice')
plt.close()

fig = plt.figure(figsize=(16, 12))
plt.plot(temps, diffslice)
plt.title('Difference from input at width=0.1', fontsize=28)
plt.xlabel('Input log(T)', fontsize=24)
plt.ylabel('Difference (%)', fontsize=24)
plt.savefig('/home/drew/Dropbox/euroscipy/diffslice')
plt.close()

fig = plt.figure(figsize=(16, 12))
ax = fig.add_subplot(1, 1, 1)
plt.plot(temps, np.log10(fitslice))
plt.title('Goodness-of-fit at width=0.1', fontsize=28)
plt.xlabel('Input log(T)', fontsize=24)
plt.ylabel('log(Goodness-of-fit)', fontsize=24)
plt.savefig('/home/drew/Dropbox/euroscipy/fitslice')
plt.close()"""

# Run synthetic data throguh 3param tempmap method

# Somehow display the results.
