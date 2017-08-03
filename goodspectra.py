### Goodman Spectrograph image reduction and spectrum extraction (with cosmic ray removal)  ###
###		Joshua Reding UNC-CH 7/27/17		###

import numpy as np
import os
import pyfits
import glob
from astropy import stats
import matplotlib.pyplot as plt
import lacosmic
from collections import Counter

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

def norayjose(img):
    print ''
    print 'Finding cosmic rays in ', img
    datalist = pyfits.open(img)
    data = np.float32(datalist[0].data)
    gain = datalist[0].header['GAIN'] #1.33 from 2017-06-07
    rdnoise = datalist[0].header['RDNOISE']
    c = lacosmic.lacosmic(data, contrast=2, cr_threshold=30, neighbor_threshold=.5, readnoise=rdnoise, effective_gain=gain, maxiter=4)
    return c[0]
###################################################
oldfiles = glob.glob('clean*.fits')
for o in oldfiles:
    os.remove(o)

objname, grating = object()
print "Finding spectra for "+objname+" taken with "+grating+"-line grating..."

biasnames = glob.glob('*Zero*.fits')
bias_med = medcomb(biasnames)[:, 50:]
flatnames = glob.glob('*QuartzFlat*.fits')
flat_med = medcomb(flatnames)[:, 50:]

imnames = glob.glob('*'+objname+'_'+grating+'.fits')
imnames = sorted(imnames, key=lambda imsa: int(imsa.split('_')[0]))
pxscl = 0.15 #"/pix
if objname == 'WD_J2350':
    sumspec = []
    n = 0
for img in imnames:
    im_clean = norayjose(img)[:, 50:]
    imdata = (im_clean - bias_med) / (flat_med - bias_med) * np.average(flat_med - bias_med)
    hdu = pyfits.PrimaryHDU(imdata)
    hdu.writeto('clean_'+img.split('.')[0]+'_red.fits', clobber=True)
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

    plt.figure()
    plt.xlabel('Pixel')
    plt.ylabel('Flux')
    plt.plot(spectrum)
    plt.show()

#     #Extract lamp spectra
#     hgarlist = glob.glob('*HgAr*.fits')
#     if hgarlist != []:
#         hgar = medcomb(hgarlist)[:, 50:]
#         #Find line locations
#         lx = np.arange(np.shape(hgar)[1])
#         max1 = 3600 + np.argmax(hgarspec[3600:3900])
#         hgar1 = np.sum(hgarspec[max1-4:max1+4]*lx[max1-4:max1+4])/np.sum(hgarspec[max1-4:max1+4])
#         max2 = 4000 + np.argmax(hgarspec[4000:])
#         hgar2 = np.sum(hgarspec[max2-4:max2+4]*lx[max2-4:max2+4])/np.sum(hgarspec[max2-4:max2+4])
#         wavelengths = [696.54, 706.72]
#         hgarlines = [hgar1, hgar2]
#     else:
#         print "No HgAr lamp calibration images found"
#     felist = glob.glob('*Fe*.fits')
#     if felist != []:
#         fe = medcomb(felist)[:, 50:]
#     else:
#         print "No Fe lamp calibration images found"
#
#     if objname == 'WD_J2350':
#         if n < 4:
#             sumspec.append(spectrum)
#         n += 1
#
# if objname == 'WD_J2350':
#     sumspec = stats.sigma_clip(np.sum(np.array(sumspec), axis=0), sigma=5)
#
#     plt.figure()
#     plt.plot(sumspec)
#     plt.show()
