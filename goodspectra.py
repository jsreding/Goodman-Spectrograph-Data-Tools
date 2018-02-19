### Goodman Spectrograph image reduction and spectrum extraction (with cosmic ray removal)  ###
### Joshua Reding UNC-CH 7/27/17   ###

import numpy as np
import os
import pyfits
import glob
from astropy import stats
import matplotlib.pyplot as plt
import lacosmic
from collections import Counter
from lmfit.models import GaussianModel

###################################################
def object():
    obj = raw_input("Prefix (e.g. WD): ")
    name = raw_input("RA (e.g. 0037): ")
    objname = obj+'_'+name
    setup = raw_input("Setup [930B, 930R, 1200Ha, 1200Ca]: ")
    return objname, setup

def medcomb(names):
    ims = []
    for n in range(len(names)):
        ims.append(pyfits.open(names[n])[0].data)
    im_med = np.median(np.array(ims), axis=0)
    return im_med

def norayjose(img, crt):
    print ''
    print 'Finding cosmic rays in ', img
    datalist = pyfits.open(img)
    hdr = datalist[0].header
    data = np.float32(datalist[0].data)
    gain = hdr['GAIN'] #1.33 from 2017-06-07
    rdnoise = hdr['RDNOISE']
    c = lacosmic.lacosmic(data, contrast=2, cr_threshold=int(crt), neighbor_threshold=.5, readnoise=rdnoise, effective_gain=gain, maxiter=4)
    return c[0][:, 50:], hdr

def arcgauss(r, p):
    p = int(p)
    r = int(r)
    l = GaussianModel(prefix='l_')
    pars = l.guess(arcdata[r, p-5:p+5], x=x[p-5:p+5])
    pars.update(l.make_params())
    pars['l_center'].set(p, min=p-5, max=p+5)
    pars['l_sigma'].set(0.5, max=np.sqrt(3))
    pars['l_amplitude'].set(50, min=0)
    mod = l
    # plt.plot(x[p-5:p+5], arcdata[r, p-5:p+5], label='data')
    out = mod.fit(arcdata[r, p-5:p+5], pars, x=x[p-5:p+5])
    comps = out.eval_components(x=x[p-5:p+5])
    # print(out.fit_report(min_correl=0.5))
    # plt.plot(x[p-5:p+5], out.best_fit, 'r-', label='best fit')
    # plt.plot(x[p-5:p+5], comps['l_'], 'b--', label='line')
    # plt.legend()
    # plt.show()
    return out.params['l_center']
###################################################
oldfiles = glob.glob('clean*.fits')+glob.glob('spec*.fits')
for o in oldfiles:
    os.remove(o)

objname, setup = object()
print "Finding spectra for "+objname+" taken in "+setup+"setup..."

biasnames = glob.glob('*Zero*.fits')
bias_med = medcomb(biasnames)[:, 50:]
flatnames = glob.glob('*Flat*.fits')
flat_med = medcomb(flatnames)[:, 50:]

imnames = glob.glob('*'+objname+'_'+grating+'.fits')
imnames = sorted(imnames, key=lambda imsa: int(imsa.split('_')[0]))
pxscl = 0.15 #"/pix
print ""
crt = raw_input("Cosmic Ray Threshold? ")

for img in imnames:
    im_clean, hdr = norayjose(img, crt)
    imdata = (im_clean - bias_med) / (flat_med - bias_med) * np.average(flat_med - bias_med)
    hdu = pyfits.PrimaryHDU(imdata)
    hdu.header = hdr
    hdu.header['BZERO'] = 0
    hdu.writeto('clean_'+img.split('.')[0]+'.fits', clobber=True)
    #Find midpoint of galaxy
    mid = np.zeros(np.shape(imdata)[1])
    for i in range(np.shape(imdata)[1]):
        mid[i] = np.argmax(imdata[:, i])
    ymid = int(Counter(mid).most_common(1)[0][0])
    #Find max Y across all X
    x = np.arange(np.shape(imdata)[1])
    y = np.zeros(np.shape(imdata)[1])
    for i in range(np.shape(imdata)[1]):
    	y[i] = ymid-5 + np.argmax(imdata[ymid-7:ymid+7, i])
    #3rd-order polynomial
    sp = np.poly1d(np.polyfit(x, y, 3))

    #Extract spectrum
    spectrum = np.zeros(np.shape(imdata)[1])
    for i in range(np.shape(imdata)[1]):
            x = np.array(np.ndarray.tolist(np.arange(20))+np.ndarray.tolist(np.arange(40,60)))
            y = np.array(np.ndarray.tolist(imdata[int(sp(i)-30):int(sp(i)-10), i]) + np.ndarray.tolist(imdata[int(sp(i)+10):int(sp(i)+30), i]))
            sk = np.poly1d(np.polyfit(x, y, 3))
            spectrum[i] = np.sum(imdata[int(sp(i)-10):int(sp(i)+10), i] - sk(np.arange(20,40)))
    spectrum = stats.sigma_clip(spectrum, sigma=5)

    #Calibrate with sky lines (Monte Carlo)

    mcpoly = np.poly1d(np.mean(mcpoly, axis=1))
    wvx = mcpoly(x)

    final = np.array([wvx, spectrum])

    plt.figure()
    plt.title(img)
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Flux (counts)')
    plt.plot(final[0], final[1])
    plt.show()
    #
    # plt.figure()
    # plt.scatter(mcl, wvs, label='sky lines')
    # plt.plot(x, wvx, label='fit')
    # plt.legend()
    # plt.show()
    # print wvs - mcpoly(mcl)

    hdu = pyfits.PrimaryHDU(final)
    hdu.header = hdr
    hdu.header['BZERO'] = 0
    hdu.writeto('spec_'+img.split('.')[0]+'.fits', clobber=True)
