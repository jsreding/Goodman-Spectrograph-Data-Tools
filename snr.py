#call as 'python snr.py diam time thrpt airm mag pxscl numpix srcltfrac QE drkcurr rdout'

import numpy as np
import sys

d = float(sys.argv[1]) #telescope diameter (m)
t = float(sys.argv[2]) #exposure time (s)
thrpt = float(sys.argv[3]) #throughput (0 - 1)
am = float(sys.argv[4]) #airmass
m = float(sys.argv[5]) #source magnitude
pxscl = float(sys.argv[6]) #pixel scale (asec/px)
npix = float(sys.argv[7]) # # of pixels
frac = float(sys.argv[8]) #fraction of source light in pixels (0 - 1)
qe = float(sys.argv[9]) #quantum efficiency (0 - 1)
dk = float(sys.argv[10]) #dark current
rd = float(sys.argv[11]) #readout (e-)

nphot = 717. # # of photons; given (photons/s/angstroms/cm^2)
bw = 1330. #bandwidth; given (angstroms)
sky = 20.9 #sky brightness (magnitudes/arcsecond^2)

mag0F = nphot * bw
ext = 0.13 * am

def area(diam):
	return np.pi * ((diam * 100.)/2)**2
def bg(diam, time, npx, pscl, thrupt, QE):
	return mag0F * 10**(-0.4 * sky) * QE * thrupt * area(diam) * (npx * pscl**2.) * time
def src(diam, time, npx, pscl, thrupt, QE, airm, frc):
	return mag0F * 10**(-0.4 * (m + ext)) * QE * thrupt * area(diam) * frc * time
def dark(drk, npx, time):
	return drk * npx * time
def rdout(rdot, npx):
	return npx * (rdot**2.)

def SNR(diam, time, npx, pscl, thrupt, QE, airm, frc, drk, rdot):
	sig = src(diam, time, npx, pscl, thrupt, QE, airm, frc)
	nois = np.sqrt(bg(diam, time, npx, pscl, thrupt, QE) + src(diam, time, npx, pscl, thrupt, QE, airm, frc) + dark(drk, npx, time) + rdout(rdot, npx))
	return np.round(sig/nois, 2)

print 'Source photons =', src(d, t, npix, pxscl, thrpt, qe, am, frac)/t
print 'Background photons =', bg(d, t, npix, pxscl, thrpt, qe)/t
print 'Photon noise =', src(d, t, npix, pxscl, thrpt, qe, am, frac)
print 'Background noise =', bg(d, t, npix, pxscl, thrpt, qe)
print 'Dark Current noise =', dark(dk, npix, t)
print 'Readout noise =', rdout(rd, npix)
print 'SNR =', SNR(d, t, npix, pxscl, thrpt, qe, am, frac, dk, rd)
