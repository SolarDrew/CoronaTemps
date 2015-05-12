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
import sunpy
from sunpy.map import Map
from temperature import TemperatureMap
from utils import gaussian, load_temp_responses
from os import path, makedirs
import subprocess32 as subp

# Define CoronaTemps home folder and output folder
CThome = path.join(path.expanduser('~'), 'CoronaTemps')
outdir = path.join(CThome, 'validation')
if not path.exists(outdir): makedirs(outdir)

# Define parameter ranges
temps = np.arange(4.6, 7.405, 0.005)
widths = np.arange(0.01, 0.605, 0.005) # Just copying Aschwanden's range here
heights = np.array([25])#, 25, 35])
print heights
n_temps = len(temps)
n_widths = len(widths)
n_heights = len(heights)
n_vals = n_temps * n_widths * n_heights
print n_temps, n_widths, n_heights, n_vals, n_vals * 6

# Create model DEMs and synthetic emission
emission = np.zeros((6, n_temps, n_heights, n_widths))
logt = np.arange(0, 15.05, 0.05)
resp = load_temp_responses()
delta_t = logt[1] - logt[0]
for t, meantemp in enumerate(temps):
    for h, height in enumerate(heights):
        for w, width in enumerate(widths):
            dem = gaussian(logt, meantemp, width, height)
            f = resp * dem
            emission[:, t, h, w] = np.sum(f, axis=1) * delta_t

emission = emission / emission[2, :, :, :]

# Load AIA response functions
resp = load_temp_responses()

# Load unnessecary map for its metadata
voidmap = Map(sunpy.AIA_171_IMAGE)
mapmeta = voidmap.meta

# Run synthetic data through 1param tempmap method
for w, wid in enumerate(heights):#widths):
    print 'Width:', wid
    fig = plt.figure(figsize=(30, 12))
    for wl, wlength  in enumerate(['94', '131', '171', '193', '211', '335']):
        emiss = Map(emission[wl, :, w, :], mapmeta)
        emiss.cmap = sunpy.cm.get_cmap('sdoaia{}'.format(wlength))
        emiss.meta['naxis1'] = emiss.shape[1]
        emiss.meta['naxis2'] = emiss.shape[0]
        emiss.meta['cdelt1'] = widths[1] - widths[0]
        emiss.meta['cdelt2'] = temps[1] - temps[0]
        emiss.meta['crval1'] = widths[0]
        emiss.meta['crval2'] = temps[0]
        emiss.meta['crpix1'] = 0.5
        emiss.meta['crpix2'] = 0.5
        fits_dir = path.join(CThome, 'data', 'synthetic', wlength)
        if not path.exists(fits_dir): makedirs(fits_dir)
        emiss.save(path.join(fits_dir, 'model.fits'))
        fig.add_subplot(1, 6, wl+1)
        emiss.plot()#vmin=0.0, vmax=1.0)
        plt.title('{}/171'.format(wlength))
        plt.xlabel('Input width')
        plt.ylabel('Input log(T)')
        plt.colorbar()
    plt.savefig(path.join(outdir, 'model_emission_h={}'.format(wid)))
    plt.close()
    
    images = [Map(emission[i, :, w, :], mapmeta) for i in range(6)]
    cmdargs = "mpiexec -n 1 python {} model {} {} {} {} {} {}".format(
        path.join(CThome, 'create_tempmap.py'), 1, CThome, None,
        None, True, False).split()
    status = subp.call(cmdargs)
    newmap = Map(path.join(CThome, 'temporary.fits'))
    subp.call(["rm", path.join(CThome, 'temporary.fits')])
    data, meta = newmap.data, newmap.meta

    truetemp = np.array(list(temps)*n_widths).reshape((n_widths, n_temps)).T
    diff = (abs(truetemp - data) / truetemp) * 100
    
    fig = plt.figure(figsize=(24, 12))
    fig.add_subplot(1, 3, 1)
    plt.imshow(data, cmap='coolwarm', origin='lower', aspect='auto',
               extent=[widths[0], widths[-1], temps[0], temps[-1]],
               interpolation='none', vmin=5.6, vmax=7.0)
    plt.colorbar()
    plt.title('Solution log(T)', fontsize=28)
    plt.ylabel('Input log(T)', fontsize=24)
    plt.xlabel('Input width', fontsize=24)
    
    fig.add_subplot(1, 3, 2)
    plt.imshow(diff, cmap='RdYlGn_r', origin='lower', aspect='auto',
               extent=[widths[0], widths[-1], temps[0], temps[-1]],
               interpolation='none')
    plt.colorbar()
    plt.title('Difference from input (%)', fontsize=28)
    plt.xlabel('Input width', fontsize=24)
    
    fig.add_subplot(1, 3, 3)
    plt.imshow(np.log10(fits), cmap='cubehelix', origin='lower', aspect='auto',
               extent=[widths[0], widths[-1], temps[0], temps[-1]],
               interpolation='none')#, vmin=fits.mean()-(2*fits.std()),
    plt.colorbar()
    plt.title('log(Goodness-of-fit)', fontsize=28)
    plt.xlabel('Input width', fontsize=24)
    plt.savefig('plots2/tempsolutions_em={:.2f}.pdf'.format(wid))
    plt.close()
    
w = np.where((widths > 0.097)*(widths < 0.103))
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
plt.close()

# Run synthetic data throguh 3param tempmap method

# Somehow display the results.
