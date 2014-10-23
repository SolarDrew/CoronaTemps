# -*- coding: utf-8 -*-
"""
Created on Thu Oct 23 13:53:36 2014

@author: drew
"""

import os
import sys
import datetime as dt
import sunpy
from sunpy.net import vso
from astropy import units as u

date = sunpy.time.parse_time(sys.argv[1])
wlen = sys.argv[2]
data_dir = sys.argv[3]

fname = os.path.join(data_dir) # + stuff

# Check if file exists
if os.path.isfile(fname):
    sys.exit()

# Download data if not found
client = vso.VSOClient()
# Wavelength value for query needs to be an astropy Quantity
wquant = u.Quantity(value=int(wlen), unit='Angstrom')
qr = client.query(vso.attrs.Time(date, date + dt.timedelta(minutes=1)),
                  vso.attrs.Wave(wquant, wquant),
                  vso.attrs.Instrument('aia'),
                  vso.attrs.Provider('JSOC'))

res = client.get(qr, path=data_dir+'{file}',site='NSO',
                 methods=['URL_FILE_Rice']).wait()