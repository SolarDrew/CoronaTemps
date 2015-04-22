import os
import glob
import subprocess
from time import sleep
from sunpy.time import parse_time as parse

wlens = ['94', '131', '171', '193', '211', '335']

def download(date):
    date = parse(date)
    # Loop through wavelengths
    for w, wlen in enumerate(wlens):
        data_dir = '/whatever/'
        fits_dir = os.path.join(data_dir, '{}/{:%Y/%m/%d}/'.format(wlen, date))
        searchname = os.path.join(fits_dir, '*{:%H:%M}*fits'.format(date))
        countfiles = glob.glob(searchname)
        if countfiles == 1:
            f = open('runfile_{}'.format(w), 'w')
            f.write("python getdata.py {} {} {}".format(date, wlen, fits_dir))
            f.close()
        # Download data if not enough found
        subprocess.call("qsub runfile_{}".format(w))

    running = True
    iterations = 0
    while running:
        if subprocess.check_output('qstat') != '':
            running = False
        elif iterations > 7200:
            
        sleep(1)
        iterations += 1
    
    return