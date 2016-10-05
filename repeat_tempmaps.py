# -*- coding: utf-8 -*-
"""
Code to make temperature maps for a given range of times and plot a histogram of
temperatures over that time range.

Created on Mon Nov 25 15:05:50 2013

@author: ajl7
"""

from matplotlib import use
use('agg')
import sunpy
from sunpy.net import hek
from sunpy.time import parse_time as parse
from sunpy.time.timerange import TimeRange as tr
import numpy as np
import datetime as dt
from sys import path
path.append('/imaps/holly/home/ajl7/CoronaTemps/')
from temperature import TemperatureMap

class DownloadError(Exception):
    def __init__(self, msg):
        self.msg = msg

def goes2flux(flareclass):
    flareclass = flareclass.upper()
    conversions = {'A': 1.0e-8, 'B': 1.0e-7, 'C': 1.0e-6, 'M': 1.0e-5,
                   'X': 1.0e-4}
    flux = float(flareclass[1:4]) * conversions[flareclass[0]]
    
    return flux


def repeat(start, end, roi_times=None, timeres=2, coords=None, ar=None,
           split_temps=None, em_wlen=None, plotminmax=False, plotstd=False,
           hist_type='plain', loaddata=False):#, output=None):
    """
    Wrapper around the main temperature map code.
    Runs create_tempmap repeatedly with the same parameters, allowing the user
    to easily create temperature maps for a particular region over an arbitrary
    period and with any temporal resolution up to that of AIA.
    
    Also calculates a histgram of the temperatures for each map and stacks them
    to create an image showing how the distribution of temperatures changes over
    the chosen time range.
    
    Parameters
    ----------
    start: str, datetime.datetime
        Start time of event being viewed. Not necessarily time of first 
        temperature map (see timerange_buffer below).

    end: str, datetime.datetime
        End time of event. Not necessarily time of last temperature map (see 
        timerange_buffer).

    timeres: int, float
        Temporal resolution of time range in hours.
    """
    #if isinstance(output, str):
    #    from matplotlib import use
    #    use(output)
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates

    #loaddata = False

    print start, end
    start, end = parse(start), parse(end)
    
    s = []
    t = []
    p = []
    
    #if flares == []:
    #    return s, t, p

    timerange = tr(start, end)
    delta = dt.timedelta(hours=timeres)
    ntimes = int(timerange.seconds()/delta.total_seconds())
    times = [time.start() for time in timerange.split(ntimes)]
    
    ntemps = 141
    tempsovertime = np.zeros((ntemps, ntimes))
    
    means = np.zeros(len(times))
    p95s = np.zeros(len(times))
    loopmeans = np.zeros(len(times))
    if plotminmax:
        maxes = np.zeros(len(times))
        mins = np.zeros(len(times))
    if plotstd:
        stds = np.zeros(len(times))
        loopstds = np.zeros(len(times))
    if em_wlen:
        meanem = np.zeros(len(times))
        if plotminmax:
            maxem = np.zeros(len(times))
            minem = np.zeros(len(times))
        if plotstd:
            stdem = np.zeros(len(times))

    for i, date in enumerate(times):
        data_only = True
        try:
            if ar == 'all':
                plotar = None
            else:
                plotar = ar
            results = output_maps(date, plotar, coords, 'data', split_temps,
                                  subimsize=50, calc_em=em_wlen, data_only=data_only)#True)#, linear=True)
            if isinstance(results, tuple):
                tempmap, emmap = results
            else:
                tempmap = results
            data = tempmap.data
        except DownloadError as de:
            data = np.zeros((512, 512))
            print de.msg
        except:
            print 'KHAAAAAAAN! Some part of the temperature-plotting process failed.'
            raise
            data = np.zeros((512, 512))
            if em_wlen:
                emmap = np.zeros((512, 512))
        
        t.append(np.nanmean(data))
        p.append(np.nanmax(data))
        
        data = data.flatten()
        data2 = data.copy()
        data2[data == 0.0] = np.NaN
        data2 = data2[np.isfinite(data)]
        data2.sort()
        temps, bins = np.histogram(data, bins=ntemps, density=False, range=(5.6, 7.0))
        temps = (temps/float(data.size))*100.0
        tempsovertime[:, i] = temps

        #loops = data[data >= split_temps]
        #data = data[data < split_temps]

        means[i] = np.nanmean(data2)
        try:
            p95s[i] = data2[round(0.95 * len(data2))-1]
        except IndexError:
            p95s[i] = np.NaN
        #loopmeans[i] = np.nanmean(loops)
        if plotminmax:
            maxes[i] = np.nanmax(data)
            mins[i] = np.nanmin(data)
            if em_wlen:
                maxem[i] = np.nanmax(emmap)
                minem[i] = np.nanmin(emmap)
        if plotstd:
            stds[i] = np.nanstd(data)
            if em_wlen:
                stdem[i] = np.nanstd(emmap)
            #loopstds[i] = np.nanstd(loops)
    
    tempsovertime[tempsovertime <= 0.1] = np.nan

    xmin, xmax = mdates.datestr2num([str(start), str(end)])
    fig = plt.figure(figsize=(36, 18))
    ax = fig.add_subplot(111, axisbg='k')
    plot_title = 'Temperature distribution of corona\n{:%Y/%m/%d %H:%M} - {:%Y/%m/%d %H:%M}'.format(start, end)
    if roi_times:
        plot_title += '\nRegion observed: {:%Y/%m/%d %H:%M} - {:%Y/%m/%d %H:%M}'.format(*roi_times)
    plt.title(plot_title)
    if hist_type == 'plain':
        plt.imshow(tempsovertime[30:106, :], extent=[xmin, xmax, 5.9, 6.65],
                   aspect='auto', interpolation='none', origin='lower',
                   cmap='coolwarm', vmin=np.nanmin(tempsovertime[65:106, :]),
                   vmax=np.nanmax(tempsovertime[65:106, :]))
    elif hist_type == 'loops':
        plt.imshow(tempsovertime[65:106, :], extent=[xmin, xmax, 6.25, 6.65],
                   aspect='auto', interpolation='none', origin='lower',
                   cmap='coolwarm', vmin=np.nanmin(tempsovertime[65:106, :]),
                   vmax=np.nanmax(tempsovertime[65:106, :]))
    elif hist_type == 'full':
        plt.imshow(tempsovertime, extent=[xmin, xmax, 5.6, 7.0],
                   aspect='auto', interpolation='none', origin='lower',
                   cmap='coolwarm', vmin=np.nanmin(tempsovertime),
                   vmax=np.nanmax(tempsovertime))
    plt.tight_layout()
    ax.xaxis_date()
    fig.autofmt_xdate()
    plt.colorbar(orientation='horizontal')
    plt.savefig('/media/huw/temp-time_hists/distribution_over_time_{}'.format(ar))
    plt.close()


    means[np.where(means == 0.0)] = np.nan
    if plotstd:
        stds[np.where(stds == 0.0)] = np.nan
        loopstds[loopstds == 0.0] = np.nan

    try:
        tnums = mdates.date2num([ti for ti in times])
        print maxes
        print len(maxes)
        fig = plt.figure(figsize=(18, 14))
        ax = fig.add_subplot(111)
        plt.title('Variation of temperature over time; AR{}'.format(ar), 
                  fontsize=32)
        plt.plot(tnums, maxes, label='Maximum temperature', color='red')
        plt.axhline(np.nanmean(maxes))
        print tnums
        print len(tnums)
        ax.xaxis_date()
        fig.autofmt_xdate()
        plt.legend(loc=4, fontsize=16)
        plt.xlabel('Date', fontsize=24)
        plt.ylabel('log(T)', fontsize=24)
        #plt.savefig('/media/huw/temp_plots/temp_plot_{}_b'.format(ar))
        plt.savefig('/home/drew/Dropbox/ARs/temps_{}_b'.format(ar))
        plt.close()

        """diff = ((maxes-p95s)/p95s)*100.0
        fig = plt.figure(figsize=(18, 14))
        ax = fig.add_subplot(1, 1, 1)
        plt.title('Percentage difference between max and 95th %-ile; AR{}'.format(ar), 
                  fontsize=32)
        plt.plot(tnums, diff, color='black')
        plt.scatter(fldates, [np.nanmax(diff)]*len(fldates))
        for flare in flares:
            ax.text(sunpy.time.parse_time(flare['event_peaktime']), np.nanmax(diff)+0.01, flare['fl_goescls'][0])
        ax.xaxis_date()
        fig.autofmt_xdate()
        plt.xlabel('Date', fontsize=24)
        plt.ylabel('log(T)', fontsize=24)
        #plt.savefig('/media/huw/temp_plots/temp_plot_{}'.format(ar))
        plt.savefig('Dropbox/ARs/diffs_{}'.format(ar))
        plt.close()"""
        
    except:# ValueError:
        print "Can't plot the temperature graph because matplotlib is being a whiney douche"
        print tnums
        raise

    return s, t, p, times
