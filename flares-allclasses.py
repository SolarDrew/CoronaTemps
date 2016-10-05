# -*- coding: utf-8 -*-
"""
Created on Tue Dec 09 18:15:45 2014

@author: Drew Leonard
"""

from matplotlib import use
use('pdf')
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cmx
import matplotlib.colors as colours
from sunpy import wcs
from sunpy.net import hek
from sunpy.time import parse_time as parse
from sunpy.time.timerange import TimeRange as tr
import datetime as dt
from os.path import join, expanduser
import sys
from sys import path
path.append(expanduser(join('~', 'CoronaTemps')))
from temperature import TemperatureMap as tmap
from astropy import units


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


def p95(data):
    return np.percentile(data, 95)


def p5(data):
    return np.percentile(data, 5)

fastdata = join('/fastdata', 'sm1ajl', 'thesis')

parameter = 'mean'
functions = {'mean': np.nanmean, 'max': np.nanmax, 'min': np.nanmin, 
             '95th %-ile': p95, '5th %-ile': p5}
savedir = join(fastdata, 'plots', 'chapter6', 'all_classes', '{}'.format(parameter.replace(' ', '_')))

start = parse('2011-02-01')
end = parse('2011-04-01')

client = hek.HEKClient()
flares = client.query(hek.attrs.Time(start, end),
                      hek.attrs.EventType('FL'))

flares = [fl for fl in flares if (fl['ar_noaanum'] > 11137 and 
                                  fl['ar_noaanum'] < 11184)]
ar_rad = 75
ar_temps_fltime = []
ar_temps_1 = []
ar_temps_10 = []
ar_temps_30 = []
fl_classes = []

flarelist = open(join(savedir, "flarelist.txt"), "w")

absfig, (axa1, axa2, axa3) = plt.subplots(3, 2, sharex='col', sharey='row', figsize=(16, 24))
axa1[0].set_title('A, B and C class flares')
axa1[1].set_title('M and X class flares')
axa1[0].set_ylabel('{} log(T)'.format(parameter.title()))
axa2[0].set_ylabel('Running difference {} log(T)'.format(parameter))
axa2[0].axhline(0, linestyle='--')
axa2[1].axhline(0, linestyle='--')
axa3[0].set_ylabel('log(T) difference from flare onset time')
axa3[0].set_xlabel('Time (minutes before flare onset)')
axa3[1].set_xlabel('Time (minutes before flare onset)')
axa3[0].axhline(0, linestyle='--', color='black')
axa3[1].axhline(0, linestyle='--', color='black')
limits1 = (1000, -1000)
limits2 = (1000, -1000)
limits3 = (1000, -1000)

"""ratfig = plt.figure(figsize=(10, 10))
axr1 = ratfig.add_subplot(2, 2, 1)
axr1.set_title("T(t)/T(t=x)")
axr1.set_ylabel("Temperature ratio")
axr2 = ratfig.add_subplot(2, 2, 2)
axr2.set_title("T(t)/T(t=-1)")
axr3 = ratfig.add_subplot(2, 2, 3)
axr3.set_title("T(t)/T(t=-10)")
axr3.set_ylabel("Temperature ratio")
axr3.set_xlabel("Time (minutes)")
axr4 = ratfig.add_subplot(2, 2, 4)
axr4.set_title("T(t)/T(t=-30)")
axr4.set_xlabel("Time (minutes)")"""

# Set up some colourmap stuff for line-plotting later
cmap = cm = plt.get_cmap('afmhot')
cNorm  = colours.Normalize(vmin=0, vmax=1)
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cmap)
flarecolours = [{'A': 0.3, 'B': 0.5, 'C': 0.7}, {'M': 0.4, 'X': 0.6}]

for flare in flares:
  try:
    flaretime = parse(flare['event_starttime'])
    starttime = flaretime-dt.timedelta(hours=0.5)
    timerange = tr(starttime, flaretime)

    region = client.query(hek.attrs.EventType('AR'),
                          hek.attrs.Time(flaretime-dt.timedelta(minutes=5), 
                                         flaretime))
    if flare['ar_noaanum'] == 11149:
        print 'HEK is fucked, using 11147'
        flare['ar_noaanum'] = 11147
    region = [r for r in region if r['ar_noaanum'] == flare['ar_noaanum']]
    if isinstance(region, list):
        try:
            region = region[0]
        except:# IndexError:
            print "An error occured for active region AR{}".format(flare['ar_noaanum'])
            continue
    
    # Define times for maps
    delta = dt.timedelta(minutes=1)
    ntimes = int(timerange.seconds()/delta.total_seconds())
    times = [time.start() for time in timerange.split(ntimes)]
    
    data_root = join(fastdata, 'data')
    maps_root = join(fastdata, 'maps')
    
    # Create empty lists to store temperature values and running differences
    # Called means because I did mean first, but used for max, percentiles,
    # etc sometimes
    means = []
    #runningdiffT = []
    for time in times:
        # Load/calculate temperature map data
        print time
        try:
            if time not in [times[0], times[-11], times[-2], times[-1]]:
                means.append(np.nan)
                continue
            data_dir = join(data_root, "{:%Y/%m/%d/}".format(time))
            maps_dir = join(maps_root, "{:%Y/%m/%d/}".format(time))
            thismap = tmap(time, data_dir=data_dir, maps_dir=maps_dir)#, verbose=True)
            thismap.save()
        
            # Crop temperature map to active region
            x, y = wcs.convert_hg_hpc(region['hgc_x'], region['hgc_y'],
                                      b0_deg=thismap.heliographic_latitude,
                                      l0_deg=thismap.carrington_longitude)
            thismap = thismap.submap([x-ar_rad, x+ar_rad], [y-ar_rad, y+ar_rad])

            # Append appropriate temperature values to list
            means.append(functions[parameter](thismap.data))
            #if len(means) > 1:
            #    runningdiffT.append(means[-1]-means[-2])
            #else:
            #    runningdiffT.append(0)
            
        except:
            #failed = True
            print "Failed", time
            means.append(np.nan)
            raise

    #print means
    # Convert time values to time before flare
    times = [(t - flaretime).total_seconds()/60 for t in times]
    # Decide line colour based on flare class
    flcl = str(flare['fl_goescls'])[0].upper()
    # Plot temperature values of AR with time for that flare
    if flcl in ['A', 'B', 'C']:
        col = 0
    else:
        col = 1
    colourVal = scalarMap.to_rgba(flarecolours[col][flcl])
    #absdiffT = [mean-means[-1] for mean in means]
    limits1 = (min(limits1[0], min(means)), max(limits1[1], max(means)))
    #limits2 = (min(limits2[0], min(runningdiffT)), max(limits2[1], max(runningdiffT)))
    #limits3 = (min(limits3[0], min(absdiffT)), max(limits3[1], max(absdiffT)))
    axa1[col].plot(times, means, color=colourVal)
    #axa2[col].plot(times, runningdiffT, color=colourVal)
    #axa3[col].plot(times, absdiffT, color=colourVal) # Absolute difference
    #ax3[col].plot(times, [((mean-means[-1])/means[-1])*100 for mean in means], color=colorVal) # Percentage difference

    """# Plot ratio of temperatures to given temperature over time
    axr1.plot(times, [mean/means[0] for mean in means], label='x = -30', color=colourVal)
    axr2.plot(times, [mean/means[-11] for mean in means], label='x = -10', color=colourVal)
    axr3.plot(times, [mean/means[-2] for mean in means], label='x = -1', color=colourVal)
    axr4.plot(times, [mean/means[-1] for mean in means], label='x = 0', color=colourVal)
    #plt.ylim(0, 1)
    #plt.legend()"""

    # Append  temperature values for final temperature map to list
    ar_temps_fltime.append(means[-1])
    ar_temps_1.append(means[-2])
    ar_temps_10.append(means[-11])
    ar_temps_30.append(means[0])
    # Append class of flare to list
    fl_classes.append(np.log10(flareclass_to_flux(str(flare['fl_goescls'])).value))

    failed = False
    flarelist.write("{} {} & {} & {} \\\\ \n".format(flaretime.date(), flaretime.time(), flare['fl_goescls'], flare['ar_noaanum']))
  except:
    print 'Failed for {} flare at {}'.format(flare['fl_goescls'], flare['event_starttime'])
    failed = True
    #raise

flarelist.close()

axa1[0].set_ylim(limits1[0]-0.02, limits1[1]+0.02)
axa1[1].set_ylim(limits1[0]-0.02, limits1[1]+0.02)
"""axa2[0].set_ylim(limits2[0]-0.005, limits2[1]+0.005)
axa2[1].set_ylim(limits2[0]-0.005, limits2[1]+0.005)
axa3[0].set_ylim(limits3[0]-0.005, limits3[1]+0.005)
axa3[1].set_ylim(limits3[0]-0.005, limits3[1]+0.005)"""
absfig.savefig(join(savedir, "allars"))#"oldallars"))
#ratfig.savefig(join(savedir, "tempratios"))
plt.close('all')

limits = (1000, -1000)
# Redefine flare colours for going on the same plot.
flarecolours = {'A': 0.2, 'B': 0.35, 'C': 0.5, 'M': 0.65, 'X': 0.8}
# Plot instantaneous temperatures of active regions for all flares against flare class
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row', figsize=(18, 18))
ax1.set_title("30 minutes before flare")
#ax1.axhline(-8, color=scalarMap.to_rgba(flarecolours['A']), linestyle='--')
ax1.axhline(-7, color=scalarMap.to_rgba(flarecolours['B']), linestyle='--')
ax1.axhline(-6, color=scalarMap.to_rgba(flarecolours['C']), linestyle='--')
ax1.axhline(-5, color=scalarMap.to_rgba(flarecolours['M']), linestyle='--')
ax1.axhline(-4, color=scalarMap.to_rgba(flarecolours['X']), linestyle='--')
ax1.scatter(ar_temps_30, fl_classes)
ax1.set_ylabel("log(flux)")
limits = (min(limits[0], min(ar_temps_30)), max(limits[1], max(ar_temps_30)))
ax2.set_title("10 minutes before flare")
#ax2.axhline(-8, color=scalarMap.to_rgba(flarecolours['A']), linestyle='--')
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
ax2.scatter(ar_temps_10, fl_classes)
limits = (min(limits[0], min(ar_temps_10)), max(limits[1], max(ar_temps_10)))
ax3.set_title("1 minute before flare")
#ax3.axhline(-8, color=scalarMap.to_rgba(flarecolours['A']), linestyle='--')
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
ax3.scatter(ar_temps_1, fl_classes)
limits = (min(limits[0], min(ar_temps_1)), max(limits[1], max(ar_temps_1)))
ax3.set_ylabel("log(flux)")
ax3.set_xlabel("{} temperature of active region".format(parameter.title()))
ax4.set_title("At time of flare")
#ax4.axhline(-8, color=scalarMap.to_rgba(flarecolours['A']), linestyle='--')
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
ax4.scatter(ar_temps_fltime, fl_classes)
limits = (min(limits[0], min(ar_temps_fltime)), max(limits[1], max(ar_temps_fltime)))
ax4.set_xlabel("{} temperature of active region".format(parameter.title()))
for axis in [ax1, ax2, ax3, ax4]:
    axis.set_xlim(limits[0]-0.02, limits[1]+0.02)
plt.savefig(join(savedir, "allflares"))
plt.close()

"""limits = (1000, -1000)
# Plot temperature differences of active regions for all flares against flare class
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, sharex='col', sharey='row', figsize=(24, 8))
ax1.set_title("T(t=30) - T(t=0)")
ax1.axhline(-8, color=scalarMap.to_rgba(flarecolours['A']), linestyle='--')
ax1.axhline(-7, color=scalarMap.to_rgba(flarecolours['B']), linestyle='--')
ax1.axhline(-6, color=scalarMap.to_rgba(flarecolours['C']), linestyle='--')
ax1.axhline(-5, color=scalarMap.to_rgba(flarecolours['M']), linestyle='--')
ax1.axhline(-4, color=scalarMap.to_rgba(flarecolours['X']), linestyle='--')
dT = [ar_temps_fltime[i]-ar_temps_30[i] for i in range(len(ar_temps_fltime))]
ax1.scatter(dT, fl_classes)
limits = (min(limits[0], min(dT)), max(limits[1], max(dT)))
ax1.set_ylabel("log(flux)")
ax1.set_xlabel("Difference from flare onset time")

ax2.set_title("T(t=10) - T(t=0)")
ax2.axhline(-8, color=scalarMap.to_rgba(flarecolours['A']), linestyle='--')
ax2.axhline(-7, color=scalarMap.to_rgba(flarecolours['B']), linestyle='--')
ax2.axhline(-6, color=scalarMap.to_rgba(flarecolours['C']), linestyle='--')
ax2.axhline(-5, color=scalarMap.to_rgba(flarecolours['M']), linestyle='--')
ax2.axhline(-4, color=scalarMap.to_rgba(flarecolours['X']), linestyle='--')
dT = [ar_temps_fltime[i]-ar_temps_10[i] for i in range(len(ar_temps_fltime))]
ax2.scatter(dT, fl_classes)
limits = (min(limits[0], min(dT)), max(limits[1], max(dT)))
ax2.set_xlabel("Difference from flare onset time")

ax3.set_title("T(t=1) - T(t=0)")
ax3.axhline(-8, color=scalarMap.to_rgba(flarecolours['A']), linestyle='--')
ax3.axhline(-7, color=scalarMap.to_rgba(flarecolours['B']), linestyle='--')
ax3.axhline(-6, color=scalarMap.to_rgba(flarecolours['C']), linestyle='--')
ax3.axhline(-5, color=scalarMap.to_rgba(flarecolours['M']), linestyle='--')
ax3.axhline(-4, color=scalarMap.to_rgba(flarecolours['X']), linestyle='--')
dT = [ar_temps_fltime[i]-ar_temps_1[i] for i in range(len(ar_temps_fltime))]
ax3.scatter(dT, fl_classes)
limits = (min(limits[0], min(dT)), max(limits[1], max(dT)))
ax3.set_xlabel("Difference from flare onset time")
for axis in [ax1, ax2, ax3]:
    axis.set_xlim(limits[0]-0.02, limits[1]+0.02)

plt.savefig(join(savedir, "allflares_diffs"))
plt.close()
"""
