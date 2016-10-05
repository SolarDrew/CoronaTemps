# -*- coding: utf-8 -*-
"""
Created on Tue Dec 09 18:15:45 2014

@author: Drew Leonard
"""

from matplotlib import use, rc
use('pdf')
rc('savefig', bbox='tight', pad_inches=0.5)
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cmx
import matplotlib.colors as colours
from matplotlib import patches
from sunpy import wcs
from sunpy.map import Map
from sunpy.net import hek
from sunpy.time import parse_time as parse
from sunpy.time.timerange import TimeRange as tr
import datetime as dt
import os
from os.path import join, exists, dirname, basename
import sys
from sys import path, argv
path.append('/home/sm1ajl/CoronaTemps/')
from temperature import TemperatureMap as tmap
from astropy import units
import glob
from scipy.io.idl import readsav


def flareclass_to_flux(flareclass):
    """
    Converts a GOES flare class into the corresponding X-ray flux.
    
    Parameters
    ----------
    flareclass : string
        The flare class to convert into X-ray flux, as a string. 
        E.g.: 'X3.2', 'm1.5', 'A9.6'.
    
    Returns
    -------
    flux : astropy.units.Quantity
        X-ray flux between 1 and 8 Angstroms as measured near Earth in W/m^2
    
    Examples
    --------
    >>> flareclass_to_flux('A1.0')
    1e-08
    >>> flareclass_to_flux('c4.7')
    4.7e-06
    >>> flareclass_to_flux('X2.4')
    0.00024

    """
    assert isinstance(flareclass, str)
    flareclass = flareclass.upper()
    conversions = {'A': 1.0e-8, 'B': 1.0e-7, 'C': 1.0e-6, 'M': 1.0e-5,
                   'X': 1.0e-4}
    fluxval = float(flareclass[1:]) * conversions[flareclass[0]]
    flux = units.Quantity(fluxval, "W/m^2")
    
    return flux


density = None#'211'

parameter = 'all'#'min'
allpars = ['min', '5th %-ile', '10th %-ile', 'mean', 'stddev', '90th %-ile', '95th %-ile', 'max',
           'n over 6.1', 'n over 6.3', 'n over 6.5']#,
#           'n 6.02-6.08']
if parameter == 'all':
    pars = allpars[:8]
else:
    pars = [parameter]

savedir = "/fastdata/sm1ajl/thesis/plots/chapter6/x_class/temperature/{}/".format(parameter.replace(' ', '_').replace('.', '_'))
print exists(savedir)
if density:
    savedir = savedir.replace('temperature', 'density_{}'.format(density))
if not exists(savedir):
    os.makedirs(savedir)

#flarefiles = glob.glob('/imaps/sspfs/archive/sdo/aia/flares/*/*.dat')
#print len(flarefiles)
start = parse('2010-07-14')
end = parse('2011-09-19')

client = hek.HEKClient()
allflares = client.query(hek.attrs.Time(start, end),
                         hek.attrs.EventType('FL'))
allflares = [fl for fl in allflares if (fl['fl_goescls'] != "" and fl['fl_goescls'][0] == 'X')]
print len(allflares)

ar_rad = 75
ar_temps_fltime = []
ar_temps_1 = []
ar_temps_10 = []
ar_temps_30 = []
fl_classes = []

flarelist = open(join(savedir, "dTdt.txt"), "w")
preamble = """\\begin{tabular}{c|c|c|c|c}
Date and time & $\\langle \\frac{\\Delta T_{5th}}{\\Delta t} \\rangle$ & $\\langle \\frac{\\Delta T_{mean}}{\\Delta t} \\rangle$ & $\\langle \\frac{\\Delta T_{95th}}{\\Delta t} \\rangle$ \\\\
\hline\n"""
flarelist.write(preamble)

tmapfig = plt.figure("tmaps", figsize=(32, 24))

figs, axa1 = [], []
for i in range(4):
    thisfig, thisax = plt.subplots(3, 2, figsize=(16, 24))
    figs.append(thisfig)
    axa1 += list(thisax.flatten())
print len(figs), len(axa1)
#absfig, axa1 = plt.subplots(4, 6, figsize=(32, 16))
#axa1 = axa1.flatten()
#axa1.set_ylabel('{} log(T)'.format(parameter.title()))
limits1 = (1000, -1000)
limits2 = (1000, -1000)
limits3 = (1000, -1000)

#absfig.suptitle(parameter)

# Set up some colourmap stuff for line-plotting later
cmap = cm = plt.get_cmap('afmhot')
cNorm  = colours.Normalize(vmin=0, vmax=1)
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cmap)
flarecolours = {'A': 0.2, 'B': 0.3, 'C': 0.4, 'M': 0.5, 'X': 0.6}

#allars_hist = np.zeros((141, 30))
dTdt = np.zeros(shape=(3, len(allflares)))#len(flarefiles)-4))

f = 0
#for file in flarefiles[:-2]:
#  flare = readsav(file)['dnow'][0]
for flare in allflares:
  #keys = flare.dtype.fields.keys()
  #keys.sort()
  #for x in keys[:12]:
  #  print x, flare[x]
  #sys.exit()
  #print '-------'
  #print flare['AR'], flare['xpos'], flare['ypos']
  #continue
  try:
    flaretime = dt.datetime(1958, 1, 1) + dt.timedelta(0, flare['starttai'])
    if flaretime.minute % 2 != 0: flaretime = flaretime.replace(minute=flaretime.minute - 1)
    starttime = flaretime - dt.timedelta(hours=0.5)
    timerange = tr(starttime, flaretime)
    x, y = flare['xpos'], flare['ypos']
    if (x, y) == (0, 0):
        continue

    # Define times for maps
    delta = dt.timedelta(minutes=2)
    ntimes = int(timerange.seconds()/delta.total_seconds())
    times = [time.start() for time in timerange.split(ntimes)]
    print times
    
    flaredir = dirname(file)
    print flaredir
    data_dir = join(flaredir, 'data/')
    maps_root = join(flaredir, 'fits/')
    paramvals_fname = join(flaredir, 'params_{}'.format(str(flaretime).replace(' ', 'T')))
    if density:
        paramvals_fname += '_density{}'.format(density)
    print 'Loading ', paramvals_fname
    
    if not exists(paramvals_fname):
      paramvals = np.zeros((11, ntimes))
      for t, time in enumerate(times):
        # Load/calculate temperature map data
        print time
        try:
            """cname = join(data_dir, '171/AIA{0:%Y%m%d}?{0:%H%M}*fits'.format(time))
            print dirname(cname), exists(dirname(cname))
            print cname, exists(cname)
            coordsmap = Map(cname)
            try:
             if isinstance(coordsmap, list): coordsmap = coordsmap[0]
            except:
             print cname, exists(cname)
             raise
            cfig, cax = plt.subplots(figsize=(16, 16))
            coordsmap.plot()
            rect = patches.Rectangle([x-ar_rad, y-ar_rad], ar_rad*2, ar_rad*2, color='blue',
                                     fill=False)
            cax.add_artist(rect)
            output_dir = join('/imaps/holly/home/ajl7/AR-tmaps/',
                              basename(flaredir), '{:%Y/%m/%d}'.format(time))
            cname = '{}/{:%Y-%m-%dT%H%M}_fulldisk'.format(output_dir, time)
            plt.savefig(cname)
            plt.close()"""

            maps_dir = join(maps_root, "{:%Y/%m/%d}/temperature/".format(time))
            thismap = tmap(time, data_dir=data_dir, maps_dir=maps_dir)#, verbose=True)#,
            #               submap=([x-ar_rad, x+ar_rad], [y-ar_rad, y+ar_rad]))
            thismap.save()
            if density:
                thismap = thismap.calculate_density(wlen=density)
            # Crop temperature map to active region
            largemap = thismap.submap([x-200, x+200], [y-200, y+200])
            thismap = thismap.submap([x-ar_rad, x+ar_rad], [y-ar_rad, y+ar_rad])
            plt.figure("tmaps")
            tmapfig.clf()
            axt1 = tmapfig.add_subplot(111, axisbg='black')
            if density:
                largemap.plot(vmin=np.nanmean(largemap.data) - (2*(np.nanstd(largemap.data))),
                              vmax=np.nanmean(largemap.data) + (2*(np.nanstd(largemap.data))))
            else:
                #largemap.data[largemap.data < 6.02] = np.NaN
                #largemap.data[largemap.data > 6.08] = np.NaN
                largemap.plot()
                plt.colorbar()
                """axt2 = tmapfig.add_subplot(122, axisbg='black')
                selmap = Map(largemap.data.copy(), largemap.meta.copy())
                selmap.cmap = largemap.cmap
                selmap.plot()
                plt.colorbar()"""
            rect = patches.Rectangle([x-ar_rad, y-ar_rad], ar_rad*2, ar_rad*2, color='white',
                                     fill=False)
            axt1.add_artist(rect)
            #axt2.add_artist(rect)
            #plt.colorbar()
            #plt.tight_layout()
            plt.title('Max: {:.2f}, Mean: {:.2f}, Min: {:.2f}'.format(thismap.max(), thismap.mean(), thismap.min()))
            output_dir = join('/imaps/holly/home/ajl7/AR-tmaps/',
                              basename(flaredir), '{:%Y/%m/%d}'.format(time))
            if density:
                output_dir = output_dir.replace('tmaps', 'nmaps_{}'.format(density))
            if not exists(output_dir):
                os.makedirs(output_dir)
            plt.savefig(join(output_dir, '{:%Y-%m-%dT%H%M}'.format(time)))
        
            data = thismap.data
            paramvals[:, t] = [np.nanmin(data),
                               np.percentile(data, 5),
                               np.percentile(data, 10),
                               np.nanmean(data),
                               np.std(data),
                               np.percentile(data, 90),
                               np.percentile(data, 95),
                               np.nanmax(data),
                               np.count_nonzero(data > 6.1),
                               np.count_nonzero(data > 6.3),
                               np.count_nonzero(data > 6.5)]
            #                   np.count_nonzero(data > 6.02 and data < 6.08)]
        except:
            print "Failed", time
            paramvals[:, t] = 0
            np.savetxt(paramvals_fname, paramvals)
            raise
      np.savetxt(paramvals_fname, paramvals)
    else:
      paramvals = np.loadtxt(paramvals_fname)

    p1 = 0
    for p, par in enumerate(pars):
        # Rename a thing so I don't have to change a load of code.
        # I'm a bad man
        parind = allpars.index(par)
        means = paramvals[parind, :]
        print parind, par
        if par == 'stddev':
            continue
        #print means
        if 'n over' not in parameter and 0 in means:
            print 'Anomalous data - skipping'
            continue
        #means = np.array([i/max(means) for i in means])
        #print means
        # Convert time values to time before flare
        times2 = [(ti - flaretime).total_seconds()/60 for ti in times]
        # Decide line colour based on flare class
        #flcl = str(flare['fl_goescls'])[0].upper()
        # Plot temperature values of AR with time for that flare
        #colourVal = scalarMap.to_rgba(flarecolours[flcl])
        #col = ['B', 'C', 'M', 'X'].index(flcl)
        """absdiffT = means[:] - means[-1]
        runningdiffT = means.copy()
        runningdiffT[0] = 0
        runningdiffT[1:] -= means[:-1]
        limits1 = (min(limits1[0], min(means)), max(limits1[1], max(means)))
        limits2 = (min(limits2[0], min(runningdiffT)), max(limits2[1], max(runningdiffT)))
        limits3 = (min(limits3[0], min(absdiffT)), max(limits3[1], max(absdiffT)))"""
        plt.sca(axa1[f])
        if par == 'mean':
            plt.errorbar(times2, means, yerr=paramvals[parind+1, :])
        else:
            plt.plot(times2, means)#, color=colourVal)
        plt.title(flaretime)
        plt.xlabel('Time before flare')
        plt.ylabel('log(T)')

        """if parameter != 'n over threshold':
            tempinds = [int(round((m - 5.6) * 100)) for m in means]
            allars_hist[tempinds, times] += 1
        print 'Histogram successful'"""
    
        # Append  temperature values for final temperature map to list
        """ar_temps_fltime.append(means[-1])
        ar_temps_1.append(means[-2])
        ar_temps_10.append(means[-11])
        ar_temps_30.append(means[0])"""
        # Append class of flare to list
        #fl_classes.append(np.log10(flareclass_to_flux(str(flare['fl_goescls'])).value))
        #print 'Lists kerfuffled'

        if 0 not in paramvals[p] and par in ['5th %-ile', 'mean', '95th %-ile']:
            maxtemps = paramvals[p]
            # Change in log(T) per minute
            dTdt[p1, f] = abs((maxtemps[1:] - maxtemps[:-1])/2.0).mean()
            p1 += 1

    values = ''.join([' & ${:e} }}$'.format(thispar) for thispar in dTdt[:, f]]).replace('e', '\\times 10^{')
    flarelist.write("{} {}{} \\\\ \n".format(flaretime.date(), flaretime.time(), values))#flare['fl_goescls'], flare['ar_noaanum']))
    #print 'Flarelist file written to'
    f += 1
  except:
    print 'Failed for {} flare at {}'.format(flare, parse(flare['starttai']))
    raise
    #continue

values = ''.join([' & ${:e} }}$'.format(thispar) for thispar in dTdt.mean(axis=1)]).replace('e', '\\times 10^{')
flarelist.write("\hline \nAvg over all regions{}\\\\ \n".format(values))
flarelist.write("\\end{tabular}")
flarelist.close()


print 'plot 1'

#axa1.set_ylim(limits1[0]-0.02, limits1[1]+0.02)
for p in axa1:
    p.set_ylim(5.6, 6.8)
for i in range(4):
    figs[i].savefig(join(savedir, "allars_{}".format(i)))
for p in axa1:
    p.set_ylim(5.9, 6.2)
for i in range(4):
    figs[i].savefig(join(savedir, "allars_close_{}".format(i)))
#absfig.savefig(join(savedir, "allars_close_{}"))
plt.close('all')

print 'plot 2'

quiet = np.load('quiet.npy')
colours = [' ', 'orange', 'darkorange', 'red']

for p, par in enumerate(['5th', 'mean', '95th']):
    fig = plt.figure(figsize=(16, 12))
    plt.bar(range(len(dTdt[p])), dTdt[p])
    meanval = quiet.mean(axis=1)[p]
    plt.axhline(meanval, linestyle='--', label='Mean', color='yellow')
    for nsigs in range(1, 4):
        plussig = meanval + (dTdt.std(axis=1)[p]*nsigs)
        label = 'Mean + ${0}\sigma$ (${1:.1%}$%)'.format(nsigs, np.count_nonzero(dTdt[p]>plussig)/float(len(dTdt[p])))
        plt.axhline(plussig, linestyle='--', label=label, color=colours[nsigs])
    plt.ylabel('$\\langle \\frac{\Delta T}{\Delta t} \\rangle$')
    plt.legend()
    plt.xlim(0, 22)
    plt.savefig('dTdt-{}'.format(par))
    plt.close()

fig = plt.figure(figsize=(32, 24))
"""plt.imshow(allars_hist[int(round((limits1[0]-5.6)*100)):int(round((limits1[1]-5.6)*100))+1, :],
           cmap='gray', extent=[-30, 0, limits1[0], limits1[1]], aspect='auto',
           interpolation='none', origin='lower')
print 'imshow successful'
plt.colorbar(orientation='horizontal')
print 'colourbar successful'
plt.title('Number of active regions at a given temperature', size=24)
print 'title successful'
plt.xlabel('Time before flare (minutes)', size=20)
print 'xlabel successful'
plt.ylabel('Temperature (log(T))', size=20)
print 'ylabel successful'
plt.tight_layout()
print 'tight_layout successful'
print savedir
print join(savedir, 'allars_hist')
plt.savefig(join(savedir, "allars_hist"))
print 'savefig successful'
plt.close()

print 'plot 3'
"""
"""limits = (1000, -1000)
# Redefine flare colours for going on the same plot.
flarecolours = {'A': 0.2, 'B': 0.35, 'C': 0.5, 'M': 0.65, 'X': 0.8}
# Plot instantaneous temperatures of active regions for all flares against flare class
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row', figsize=(18, 18))
ax1.set_title("30 minutes before flare", size=24)
ax1.axhline(-8, color=scalarMap.to_rgba(flarecolours['A']), linestyle='--')
ax1.axhline(-7, color=scalarMap.to_rgba(flarecolours['B']), linestyle='--')
ax1.axhline(-6, color=scalarMap.to_rgba(flarecolours['C']), linestyle='--')
ax1.axhline(-5, color=scalarMap.to_rgba(flarecolours['M']), linestyle='--')
ax1.axhline(-4, color=scalarMap.to_rgba(flarecolours['X']), linestyle='--')
ax1.scatter(ar_temps_30, fl_classes)#, color='green')
ax1.set_ylabel("log(flux)", size=20)
#limits = (min(limits[0], min(ar_temps_30)), max(limits[1], max(ar_temps_30)))
ax2.set_title("10 minutes before flare", size=24)
ax2.axhline(-8, color=scalarMap.to_rgba(flarecolours['A']), linestyle='--')
ax2.axhline(-7, color=scalarMap.to_rgba(flarecolours['B']), linestyle='--')
ax2.axhline(-6, color=scalarMap.to_rgba(flarecolours['C']), linestyle='--')
ax2.axhline(-5, color=scalarMap.to_rgba(flarecolours['M']), linestyle='--')
ax2.axhline(-4, color=scalarMap.to_rgba(flarecolours['X']), linestyle='--')
arcols = ['blue', 'red']
for i in range(len(ar_temps_30)):
    try:
        dx = ar_temps_10[i] - ar_temps_30[i]
        ax2.arrow(ar_temps_30[i], fl_classes[i], dx, 0,
                  length_includes_head=True, head_length=0.001,
                  head_width=0.1, color=arcols[dx > 0])
    except IndexError:
        pass
ax2.scatter(ar_temps_10, fl_classes)#, color='green')
#limits = (min(limits[0], min(ar_temps_10)), max(limits[1], max(ar_temps_10)))
ax3.set_title("1 minute before flare", size=24)
ax3.axhline(-8, color=scalarMap.to_rgba(flarecolours['A']), linestyle='--')
ax3.axhline(-7, color=scalarMap.to_rgba(flarecolours['B']), linestyle='--')
ax3.axhline(-6, color=scalarMap.to_rgba(flarecolours['C']), linestyle='--')
ax3.axhline(-5, color=scalarMap.to_rgba(flarecolours['M']), linestyle='--')
ax3.axhline(-4, color=scalarMap.to_rgba(flarecolours['X']), linestyle='--')
for i in range(len(ar_temps_10)):
    try:
        dx = ar_temps_1[i] - ar_temps_10[i]
        ax3.arrow(ar_temps_10[i], fl_classes[i], dx, 0,
                  length_includes_head=True, head_length=0.001,
                  head_width=0.1, color=arcols[dx > 0])
    except IndexError:
       	pass
ax3.scatter(ar_temps_1, fl_classes)#, color='green')
#limits = (min(limits[0], min(ar_temps_1)), max(limits[1], max(ar_temps_1)))
ax3.set_ylabel("log(flux)", size=20)
ax3.set_xlabel("{} temperature of active region".format(parameter.title()), size=20)
ax4.set_title("At time of flare", size=24)
ax4.axhline(-8, color=scalarMap.to_rgba(flarecolours['A']), linestyle='--')
ax4.axhline(-7, color=scalarMap.to_rgba(flarecolours['B']), linestyle='--')
ax4.axhline(-6, color=scalarMap.to_rgba(flarecolours['C']), linestyle='--')
ax4.axhline(-5, color=scalarMap.to_rgba(flarecolours['M']), linestyle='--')
ax4.axhline(-4, color=scalarMap.to_rgba(flarecolours['X']), linestyle='--')
for i in range(len(ar_temps_1)):
    try:
        dx = ar_temps_fltime[i] - ar_temps_1[i]
        ax4.arrow(ar_temps_1[i], fl_classes[i], dx, 0,
                  length_includes_head=True, head_length=0.001,
                  head_width=0.1, color=arcols[dx > 0])
    except IndexError:
       	pass
ax4.scatter(ar_temps_fltime, fl_classes)#, color='green')
#limits = (min(limits[0], min(ar_temps_fltime)), max(limits[1], max(ar_temps_fltime)))
ax4.set_xlabel("{} temperature of active region".format(parameter.title()), size=20)
#for axis in [ax1, ax2, ax3, ax4]:
#    axis.set_xlim(limits[0]-0.02, limits[1]+0.02)
plt.savefig(join(savedir, "allflares"))
plt.close()

print 'plot 4'

limits = (1000, -1000)
# Plot temperature differences of active regions for all flares against flare class
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, sharex='col', sharey='row', figsize=(24, 8))
ax1.set_title("T(t=30) - T(t=0)", size=24)
ax1.axhline(-8, color=scalarMap.to_rgba(flarecolours['A']), linestyle='--')
ax1.axhline(-7, color=scalarMap.to_rgba(flarecolours['B']), linestyle='--')
ax1.axhline(-6, color=scalarMap.to_rgba(flarecolours['C']), linestyle='--')
ax1.axhline(-5, color=scalarMap.to_rgba(flarecolours['M']), linestyle='--')
ax1.axhline(-4, color=scalarMap.to_rgba(flarecolours['X']), linestyle='--')
dT = [ar_temps_fltime[i]-ar_temps_30[i] for i in range(len(ar_temps_fltime))]
ax1.scatter(dT, fl_classes)#, color='green')
limits = (min(limits[0], min(dT)), max(limits[1], max(dT)))
ax1.set_ylabel("log(flux)", size=20)
ax1.set_xlabel("Difference from flare onset time", size=20)

print 'plot 5'

ax2.set_title("T(t=10) - T(t=0)", size=24)
ax2.axhline(-8, color=scalarMap.to_rgba(flarecolours['A']), linestyle='--')
ax2.axhline(-7, color=scalarMap.to_rgba(flarecolours['B']), linestyle='--')
ax2.axhline(-6, color=scalarMap.to_rgba(flarecolours['C']), linestyle='--')
ax2.axhline(-5, color=scalarMap.to_rgba(flarecolours['M']), linestyle='--')
ax2.axhline(-4, color=scalarMap.to_rgba(flarecolours['X']), linestyle='--')
dT = [ar_temps_fltime[i]-ar_temps_10[i] for i in range(len(ar_temps_fltime))]
ax2.scatter(dT, fl_classes)#, color='green')
limits = (min(limits[0], min(dT)), max(limits[1], max(dT)))
ax2.set_xlabel("Difference from flare onset time", size=20)

print 'plot 6'

ax3.set_title("T(t=1) - T(t=0)", size=24)
ax3.axhline(-8, color=scalarMap.to_rgba(flarecolours['A']), linestyle='--')
ax3.axhline(-7, color=scalarMap.to_rgba(flarecolours['B']), linestyle='--')
ax3.axhline(-6, color=scalarMap.to_rgba(flarecolours['C']), linestyle='--')
ax3.axhline(-5, color=scalarMap.to_rgba(flarecolours['M']), linestyle='--')
ax3.axhline(-4, color=scalarMap.to_rgba(flarecolours['X']), linestyle='--')
dT = [ar_temps_fltime[i]-ar_temps_1[i] for i in range(len(ar_temps_fltime))]
ax3.scatter(dT, fl_classes)#, color='green')
limits = (min(limits[0], min(dT)), max(limits[1], max(dT)))
ax3.set_xlabel("Difference from flare onset time", size=20)
for axis in [ax1, ax2, ax3]:
    axis.set_xlim(limits[0]-0.02, limits[1]+0.02)
"""
plt.savefig(join(savedir, "allflares_diffs"))
plt.close()

print 'end'
