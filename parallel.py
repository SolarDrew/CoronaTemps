from os import path
from sunpy.time import parse_time as parse

wlens = ['94', '131', '171', '193', '211', '335']

def download(date):
    date = parse(date)
    # Loop through wavelengths
    for wlen in wlens[r:r+2]:
        fits_dir = data_dir + '{}/{:%Y/%m/%d}/'.format(wlen, date)
        filename = fits_dir + 'aia*{0}*{1:%Y?%m?%d}?{1:%H?%M}*lev1?fits'.format(wlen, date)
        # Check if file exists
        if path.isfile(filename):
            # Need to make sure this doesn't cause hangups
            continue
        # Download data if not enough found
        client = vso.VSOClient()
        if temp_im == []:
            # Wavelength value for query needs to be an astropy Quantity
            wquant = u.Quantity(value=int(wlen), unit='Angstrom')
            qr = client.query(vso.attrs.Time(date,# - dt.timedelta(seconds=6),
                                             date + dt.timedelta(seconds=12)),#6)),
                              vso.attrs.Wave(wquant, wquant),
                              vso.attrs.Instrument('aia'),
                              vso.attrs.Provider('JSOC'))
            res = client.get(qr, path=fits_dir+'{file}',site='NSO',
                             methods=['URL_FILE_Rice']).wait()


if __name__ == "__main__":
    print "Hello from process ", rank, "of ", size
