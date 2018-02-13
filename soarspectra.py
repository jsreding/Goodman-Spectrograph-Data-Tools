###    Program to plot spectra from SOAR/Goodman    ###
###    Joshua Reding UNC-CH 07/11/17    ###

import pyfits
import numpy as np
import matplotlib.pyplot as plt
import sys
import glob

# Will read all .fits files if no single one specified
try:
	sys.argv[1]
except IndexError:
	files = glob.glob('*.fits')
else:
	files = [sys.argv[1]]
print files

for f in files:
    img = pyfits.open(f)
    imdata = img[0].data

    tele = img[0].header['TELESCOP']
    inst = img[0].header['INSTRUME']
    conf = img[0].header['INSTCONF']
    obj = img[0].header['OBJECT']
    ra = img[0].header['OBSRA']
    dec = img[0].header['OBSDEC']
    exp = img[0].header['EXPTIME']
    date = img[0].header['DATE-OBS'].split('T')[0]
    time = img[0].header['DATE-OBS'].split('T')[1]

    print "Filename: ", f
    print "Telescope: ", tele
    print "Instrument: ", inst
    print "Object: ", obj
    print "RA: ", ra
    print "Dec: ", dec
    print "Exposure time: ", exp
    print "Date: ", date
    print "Time: ", time
    print ""

    plt.figure()
    plt.title(obj)
    plt.xlabel('Wavelength (ang)')
    plt.ylabel('Flux')
    plt.plot(imdata[91,:], imdata[0][1])
    plt.show()
