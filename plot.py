# Python plotting script; developed by Julia Bodensteiner & extended Tomer Shenar
# Type --legend for legend
# Type --scatterr for scatter plot + errors
# Other options exist, see script below...



import matplotlib.pyplot as plt
import sys
import numpy as np
import astropy.io.fits as fits
from scipy import interpolate
import sys
from astropy.table import Table
import pandas as pd

def read_file(infile, col0=0, col1=1):
    ext = str(infile.split('.')[-1])

    # Check type of input file (fits or ascii) to read in data correctly
    if (ext == 'fits') or (ext ==  'fit'):
        wave, flux = read_fits(infile)

    elif (ext == 'gz'):
        wave, flux = read_tlusty(infile)

    elif (ext == 'dat' or ext == 'ascii' or ext == 'txt' or ext == 'nspec'):
        wave, flux = read_ascii(infile)

    elif (ext == 'tfits'):
        wave, flux = read_uvespop(infile)

    elif (ext == 'hfits'):
        wave, flux = read_hermes_normalized(infile)

    else:
        wave, flux = read_ascii(infile, col0, col1)

    return wave, flux


def read_fits(infile):
    print("%s: Input file is a fits file." % infile)

    header = fits.getheader(infile)

    if 'HIERARCH SPECTRUM EXTRACTION' in header:
        wave, flux = read_psfSpec(infile)

    elif 'INSTRUME' in header:
        ins = header['INSTRUME']
        if (ins == 'MUSE'):
            wave, flux = read_pampelMUSE(infile)

        elif (ins == 'HERMES'):
            wave, flux = read_HERMES(infile)

        elif (ins == 'FEROS'):
            wave, flux = read_FEROS(infile)
        elif (ins == 'XSHOOTER'):
            wave, flux = read_XSHOOTER(infile)
        elif (ins == 'UVES'):
            wave, flux = read_UVES(infile)
        elif (ins == 'UVES_STITCH'):
            wave, flux = read_UVES_STITCH(infile)                    
        elif (ins == 'GIRAFFE' and 'nLR' in infile):
            wave, flux = read_GIRAFFE(infile)     
        elif (ins == 'GIRAFFE'):
            #wave, flux = read_GIRAFFE(infile)               
            wave, flux = read_GIRAFFE2(infile)               
        elif (ins == 'ESPCOUDE'):
            wave, flux = read_NLA(infile)   
        elif (ins == 'COS'):
            wave, flux = read_COS(infile)               
        elif (ins == 'STIS'):
            wave, flux = read_STIS(infile)              
        else:
            print('File type unkown, trying HERMES')
            wave, flux = read_HERMES(infile)                        
    else:
        wave, flux = read_HERMES(infile)

    return wave, flux

def read_NLA(infile):
    hdul = fits.open(infile)
    header = fits.getheader(infile)
    data = hdul[0].data
    flux = data[0, :, :][0]
    raw = data[1, :, :][0]
    bkg = data[2, :, :][0]
    sigma = data[3, :, :][0]
    wl0 = header['CRVAL1']
    delt = header['CD1_1']
    pix = header['CRPIX1']
    wave = wl0 - (delt * pix - delt) + np.arange(flux.shape[0]) * delt
    return wave, flux

def read_psfSpec(infile):
    print("%s: Input file is a PSF extracted file." % infile)
    header = fits.getheader(infile)
    flux = fits.getdata(infile, ext=0)
    err = fits.getdata(infile, ext=1)
    wl0 = header['CRVAL1']  # Starting wl at CRPIX1
    delt = header['CDELT1']  # Stepwidth of wl
    pix = header['CRPIX1']  # Reference Pixel
    wave = wl0 - (delt * pix - delt) + np.arange(flux.shape[0]) * delt
    flux_out = flux, err

    return wave, flux_out


def read_pampelMUSE(infile):
    print("%s: Input file is a pampelMUSE file." % infile)
    header = fits.getheader(infile)
    flux = fits.getdata(infile, ext=0)
    #print flux
    ##print header['CRVAL1'], header['CDELT1'], header['CRPIX1']
    #err = fits.getdata(infile, ext=1)
    wl0 = header['CRVAL1']  # Starting wl at CRPIX1
    delt = header['CDELT1']  # Stepwidth of wl
    pix = header['CRPIX1']  # Reference Pixel
    #print 'hi'
    #print flux.shape[0]
    #ddasfasf
    wave = wl0 - (delt * pix - delt) + np.arange(flux.shape[0]) * delt
    #print wave
    #dsadsf
    # because pampelmuse gives the wl in m --> convert to A
    #wave = wave 
    #flux_out = flux, err

    return wave, flux




def read_MUSE(infile):
    print("%s: Input file is a MUSE file." % infile)
    header = fits.getheader(infile)
    data = fits.getdata(infile)
    wl0 = header['CRVAL1']  # Starting wl at CRPIX1
    delt = header['CDELT1']  # Stepwidth of wl
    pix = header['CRPIX1']  # Reference Pixel
    wave = wl0 - (delt * pix - delt) + np.arange(data.shape[0]) * delt

    flux = data

    return wave, flux


def read_GIRAFFE(infile):
    print("%s: Input file is a GIRAFFE file." % infile)
    header = fits.getheader(infile)
    data = fits.getdata(infile)
    wl0 = header['CRVAL1']  # Starting wl at CRPIX1
    delt = header['CDELT1']  # Stepwidth of wl
    pix = header['CRPIX1']  # Reference Pixel
    wave = wl0 - (delt * pix - delt) + np.arange(data.shape[0]) * delt
    table = Table.read(infile, hdu=1)
    flux = table['NORM_SKY_SUB_CR']
    return wave, flux

def read_GIRAFFE2(infile):
    print("%s: Input file is a GIRAFFE2 file." % infile)
    hdunum=1
    header = fits.getheader(infile)
    hdul = fits.open(infile)
    table = Table.read(infile, hdu=hdunum)
    wave = table['WAVE'][0]*10.
    flux = table['FLUX_REDUCED'][0]    
    return wave, flux



def read_HERMES(infile):
    print("%s: Input file is a HERMES file." % infile)
    header = fits.getheader(infile)
    # for files with standard wavelegth array
    if ((header['CTYPE1'] == 'WAVELENGTH') or (header['CTYPE1'] == 'AWAV')):
        flux = fits.getdata(infile)
        crval = header['CRVAL1']
        cdelt = header['CDELT1']
        naxis1 = header['NAXIS1']
        wave = crval + np.arange(0, naxis1) * cdelt

    # for files that are given in logarithmic wl array
    if (header['CTYPE1'] == 'log(wavelength)'):
        flux = fits.getdata(infile)
        crval = header['CRVAL1']
        cdelt = header['CDELT1']
        naxis1 = header['NAXIS1']
        wave = np.exp(crval + np.arange(0, naxis1)*cdelt)
    else:
        flux = fits.getdata(infile)
        crval = header['CRVAL1']
        cdelt = header['CDELT1']
        naxis1 = header['NAXIS1']
        wave = crval + np.arange(0, naxis1) * cdelt                
    return wave, flux


def read_FEROS(infile):
    print("%s: Input file is a FEROS file." % infile)
    header = fits.getheader(infile)
    flux = fits.getdata(infile)
    crval = header['CRVAL1']
    crpix = header['CRPIX1']
    cdelt = header['CDELT1']

    wave = crval - (cdelt * crpix - cdelt) + np.arange(flux.shape[0]) * cdelt
    return wave, flux


def read_XSHOOTER(infile):
    print("%s: Input file is a XSHOOTER file." % infile)
    header = fits.getheader(infile)
    hdul = fits.open(infile)
    table = Table.read(infile, hdu=1)
    wave = table['WAVE'][0]
    #flux = table['FLUX_REDUCED'][0]
    flux = table['FLUX'][0]    
    #else:
        #flux = fits.getdata(infile)
        #crval = header['CRVAL1']
        #cdelt = header['CDELT1']
        #naxis1 = header['NAXIS1']
        #wave = crval + np.arange(0, naxis1) * cdelt
    return wave, flux


def read_COS(infile):
    print("%s: Input file is a COS file." % infile)
    header = fits.getheader(infile)
    hdul = fits.open(infile)
    table = Table.read(infile, hdu=1)
    waves = table['WAVELENGTH']
    fluxes = table['FLUX'] 
    wave = []
    flux = []
    for jj in np.arange(len(waves)-1, 0, -1):
        #print jj
        for w in waves[jj]:
            wave.append(w)
        for f in fluxes[jj]:
            flux.append(f)  
    #flux = table['FLUX_REDUCED'][0]
    #flux = table['FLUX'] 
    #print wave
    #print flux
    #dsadd
    #else:
        #flux = fits.getdata(infile)
        #crval = header['CRVAL1']
        #cdelt = header['CDELT1']
        #naxis1 = header['NAXIS1']
        #wave = crval + np.arange(0, naxis1) * cdelt
    return wave, flux

def read_STIS(infile):
    print("%s: Input file is a STIS file." % infile)
    header = fits.getheader(infile)
    hdul = fits.open(infile)
    table = Table.read(infile, hdu=1)
    #print table
    wave = table['WAVELENGTH'][0]
    flux = table['FLUX'][0]
    #print waves
    #print fluxes
    #wave = []
    #flux = []
    #print 'hello'
    #print np.arange(len(waves[0])-1, 0, -1)
    #for jj in np.arange(len(waves[0])-1, 0, -1):        
        #print jj, waves[0][jj], fluxes[jj]
        #for w in waves[jj]:
            #wave.append(w)
        #for f in fluxes[jj]:
            #flux.append(f)  
    #flux = table['FLUX_REDUCED'][0]
    #flux = table['FLUX'] 
    #print wave
    #print flux
    #dsadd
    #else:
        #flux = fits.getdata(infile)
        #crval = header['CRVAL1']
        #cdelt = header['CDELT1']
        #naxis1 = header['NAXIS1']
        #wave = crval + np.arange(0, naxis1) * cdelt
    return wave, flux


def read_UVES(infile):
    print("%s: Input file is a UVES or GIRAFFE file." % infile)
    header = fits.getheader(infile)
    if 'SPEC_COM' in header:
        print('a')
        table = Table.read(infile, hdu=1)
        wave = table['wave']
        flux = table['FLUX_REDUCED']
    elif 'CRVAL1' in header:
        print('b')
        flux = fits.getdata(infile)
        crval = header['CRVAL1']
        cdelt = header['CDELT1']
        naxis1 = header['NAXIS1']
        wave = crval + np.arange(0, naxis1) * cdelt
    else:
        print('c')
        table = Table.read(infile, hdu=0)
        wave = table['WAVE'][0]
        flux = table['FLUX'][0]
    if wave[0] < 1000:  # wavelength is given in nm
        wave = wave * 10 # convert wavelength to A
    return wave, flux

def read_UVES_STITCH(infile):
    print("%s: Input file is a stitched UVES file" % infile)
    table = Table.read(infile, hdu=1)
    wave = table['WAVE']
    flux = table['FLUX']
    return wave, flux


def read_tlusty(infile):
    fill_val = 'extrapolate'
    print("%s: Input file is a TLUSTY model." % infile)

    # first the model
    wave, flux = np.loadtxt(infile, unpack=True)

    # wave array not evenly spaced => interpolate it
    s = interpolate.interp1d(wave, flux, 2, fill_value=fill_val)
    flux = s(wave)

    # now the continuum
    contfile = infile.split('.7.')[0] + '.17.gz'
    wave_cont, flux_cont = np.loadtxt(contfile, unpack=True)

    # wave array not evenly spaced => interpolate it
    s = interpolate.interp1d(wave_cont, flux_cont, 1, fill_value=fill_val)
    flux_cont = s(wave)

    # normalize the model
    flux = flux / flux_cont

    return wave, flux


def read_uvespop(infile):
    print("%s: Input file is a UVESPOP file." % infile)
    data = fits.getdata(infile)
    wave = data.field(0)
    flux = data.field(1)

    return wave, flux


def read_hermes_normalized(infile):
    print("%s: Input file is a normalized HERMES file." % infile)
    data = fits.getdata(infile)
    wave = data.field(0)
    # flux = data.field(1)
    norm = data.field(2)
    # cont = data.field(3)

    return wave, norm


def read_ascii(infile, col0=0, col1=1):
    # any type of ascii file (typically I call them .dat)
    # assumes that first column is wave and second column is flux
    print("%s: Input file is an ascii file." % infile)
    #wave, flux = np.loadtxt(infile, usecols=(0, 1), unpack=True)
    bla = (pd.read_csv(infile, header=None, delim_whitespace=True, comment='#')).values
    #bla = np.loadtxt(infile)
    ##print bla
    wave = bla[:,col0]
    flux = bla[:,col1]
    return wave, flux


def cut_muse_specrange(wave, flux):
    # cut the spectrum to MUSE wavelength
    spec_range = ((wave > 4600) & (wave < 9350))
    wave = wave[spec_range]
    flux = flux[spec_range]
    return wave, flux


def mask_muse_laser(wave, flux):
    index = [j for j in range(len(wave)) if 5746.0 < wave[j] < 6017.]
    flux[index] = np.nan

    return wave, flux


def write_pampelmuse(infile, flux_cont, err_cont, outfile):
    hdul_infile = fits.open(infile)
    hdul_new = fits.HDUList()
    primheader = hdul_infile[0].header.copy()
    secondheader = hdul_infile[1].header.copy()

    hdul_new.append(fits.PrimaryHDU(data=flux_cont, header=primheader))
    hdul_new.append(fits.ImageHDU(data=err_cont, header=secondheader))
    hdul_new.writeto(outfile)

    print("Data written to %s" % outfile)


def write_hermes(infile, flux_cont, outfile):
    hdul_infile = fits.open(infile)
    hdul_new = fits.HDUList()
    primheader = hdul_infile[0].header.copy()

    hdul_new.append(fits.PrimaryHDU(data=flux_cont, header=primheader))
    hdul_new.writeto(outfile)

    print("Data written to %s" % outfile)


def write_extracted_spectrum(header, flux, flux_err, outfilename):

    hdul_new = fits.HDUList()
    hdul_new.append(fits.PrimaryHDU(data=flux, header=header))
    hdul_new.append(fits.ImageHDU(data=flux_err))
    hdul_new.writeto(outfilename)

    print("Data written to %s" % outfilename)


def write_2Dimage(header, image, outfilename):
    hdul_new = fits.HDUList()
    hdul_new.append(fits.PrimaryHDU(data=image, header=header))
    hdul_new.writeto(outfilename)



fig, ax = plt.subplots()

Legend = False
Binning = False
ScatterErr = False
Skip = False
SkipTwice = False
Norm=False
for i in range(len(sys.argv)-1):
    print(i)
    if Skip:
        if SkipTwice:
            Skip=True
            SkipTwice=False
        else:
            Skip=False
        continue
    j = i+1
    if len(sys.argv) >1:
        if sys.argv[j] == '--legend':
            Legend = True   
            continue
        if sys.argv[j] == '--norm':
            Norm = True   
            continue    
        if sys.argv[j] == '--scatterr':
            ScatterErr = True   
            continue     
        if sys.argv[j] == '--cols':
            col0 = float(sys.argv[j+1])
            col1 = float(sys.argv[j+2])         
            Skip = True
            SkipTwice = True
            continue       
    infile = sys.argv[j]
    try:       
        try:
            wave_in, flux_in = read_file(infile, col0, col1)
        except:
            wave_in, flux_in = read_file(infile)
        flux_in = np.nan_to_num(flux_in, 1.)
    except:
        continue
    try:        
        wave_in = wave_in.astype(float)
        flux_in = flux_in.astype(float)
    except:
        wave_in, flux_in = np.loadtxt(infile, unpack=True)
    if len(flux_in) == 2:
        flux = flux_in[0]
        err = flux_in[1]
    else:
        flux = flux_in
    if Norm:
        flux /= np.mean(flux)
    # Do the plotting    
    name = str(infile).split('.fits')[0]
    if ScatterErr:
        ax.errorbar(wave_in, flux, yerr=np.loadtxt(infile)[:,2], fmt='o', linewidth=1.0, alpha=0.8, label=name)
    else:
        ax.plot(wave_in, flux, linewidth=1.0, alpha=0.8, label=name)    
if Legend:
    ax.legend()
ax.set_xlabel("Wavelength [A]")
ax.set_ylabel("Flux")
plt.show()
#np.savetxt(infile + '.txt', np.c_[wave_in, flux])
