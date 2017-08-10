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
    grating = raw_input("Grating (e.g. 1200): ")
    return objname, grating

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
###################################################
oldfiles = glob.glob('clean*.fits')+glob.glob('spec*.fits')
for o in oldfiles:
    os.remove(o)

objname, grating = object()
print "Finding spectra for "+objname+" taken with "+grating+"-line grating..."

biasnames = glob.glob('*Zero*.fits')
bias_med = medcomb(biasnames)[:, 50:]
flatnames = glob.glob('*Flat*.fits')
flat_med = medcomb(flatnames)[:, 50:]

imnames = glob.glob('*'+objname+'_'+grating+'.fits')
imnames = sorted(imnames, key=lambda imsa: int(imsa.split('_')[0]))
pxscl = 0.15 #"/pix
print ""
crt = raw_input("Cosmic Ray Threshold? ")

if objname == 'WD_J2350':
    sumspec = []
    n = 0
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
    wvs = np.array([5889.953, 5895.923, 6257.961, 6287.434, 6329.840, 6498.729, 6533.044, 6553.617, 6577.285, 6863.955, 6923.220, 6949.045, 6978.414])
    pxl = np.array([232, 251, 1399, 1493, 1628, 2170, 2280, 2347, 2423, 3358, 3554, 3640, 3736])
    N = 1000
    x = np.arange(np.shape(imdata)[1])
    skl = np.zeros(len(wvs))
    mc = np.zeros(len(wvs)).astype(np.int64)
    mcl = np.zeros(len(wvs))
    mcpoly = np.zeros((3, N))
    for n in range(N):
        for i in range(len(wvs)):
            skl[i] = pxl[i]-10 + np.mean([np.argmax(imdata[int(ymid-60),pxl[i]-5:pxl[i]+5]), np.argmax(imdata[int(ymid-40),pxl[i]-5:pxl[i]+5]), np.argmax(imdata[int(ymid-20),pxl[i]-5:pxl[i]+5]), np.argmax(imdata[int(ymid+20),pxl[i]-5:pxl[i]+5]), np.argmax(imdata[int(ymid+40),pxl[i]-5:pxl[i]+5]), np.argmax(imdata[int(ymid+60),pxl[i]-5:pxl[i]+5])])
            mc[i] = int(np.random.normal(skl[i], 4))
            mcl[i] = np.sum(spectrum[mc[i]-4:mc[i]+4]*x[mc[i]-4:mc[i]+4])/np.sum(spectrum[mc[i]-4:mc[i]+4])
        mcp = np.polyfit(mcl, wvs, 2)
        mcpoly[0, n] = mcp[0]
        mcpoly[1, n] = mcp[1]
        mcpoly[2, n] = mcp[2]
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
    # plt.scatter(mcl, wvs)
    # plt.plot(x, wvx, label='fit')
    # plt.show()
    # print wvs[goodindex] - mcfixed(mcl[goodindex])

    hdu = pyfits.PrimaryHDU(final)
    hdu.header = hdr
    hdu.header['BZERO'] = 0
    hdu.writeto('spec_'+img.split('.')[0]+'.fits', clobber=True)
