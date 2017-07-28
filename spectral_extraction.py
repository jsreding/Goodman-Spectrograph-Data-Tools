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
    gain = 1.33 #datalist[0].header['GAIN'] #1.33 from 2017-06-07
    rdnoise = datalist[0].header['RDNOISE']
    c = lacosmic.lacosmic(data, contrast=2, cr_threshold=40, neighbor_threshold=1, readnoise=rdnoise, effective_gain=gain, maxiter=4)
    return c[0]
###################################################
oldfiles = glob.glob('clean*.fits')
for o in oldfiles:
    os.remove(o)

objname = glob.glob('*WD*.fits')[0].split('_')[1]+'_'+glob.glob('*WD*.fits')[0].split('_')[2]
grating = glob.glob('*WD*.fits')[0].split('_')[3]

biasnames = glob.glob('*Zero*.fits')
bias_med = medcomb(biasnames)
flatnames = glob.glob('*QuartzFlat*.fits')
flat_med = medcomb(flatnames)

imnames = glob.glob('*'+objname+'_'+grating+'.fits')
pxscl = 0.15 #"/pix
for img in imnames:
    im_clean = norayjose(img)
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
    	y[i] = ymid-5 + np.argmax(imdata[ymid-5:ymid+5, i])
    #3rd-order polynomial
    p3 = np.poly1d(np.polyfit(x, y, 3))

    #Extract spectrum
    spectrum = np.zeros(np.shape(imdata)[1])
    lt = np.zeros(np.shape(imdata)[1])
    sk = np.zeros(np.shape(imdata)[1])
    for i in range(np.shape(imdata)[1]):
    	lt[i] = np.sum(imdata[int(p3(i)-5):int(p3(i)+5), i])
    	sk[i] = np.sum(imdata[int(p3(i)-30):int(p3(i)-20), i])
    ltc = stats.sigma_clip(lt, sigma=3)
    skc = stats.sigma_clip(sk, sigma=3)
    spectrum = ltc - skc

    plt.figure()
    plt.plot(spectrum[50:])
    plt.show()

    # #Extract lamp spectra
    # hgarlist = glob.glob('*HgAr*.fits')
    # if hgarlist != []:
    #     hgar = medcomb(hgarlist)
    #     #Find line locations
    #     lx = np.arange(np.shape(hgar)[1])
    #     max1 = 3600 + np.argmax(hgarspec[3600:3900])
    #     line1 = np.sum(hgarspec[max1-4:max1+4]*lx[max1-4:max1+4])/np.sum(hgarspec[max1-4:max1+4])
    #     max2 = 4000 + np.argmax(hgarspec[4000:])
    #     line2 = np.sum(hgarspec[max2-4:max2+4]*lx[max2-4:max2+4])/np.sum(hgarspec[max2-4:max2+4])
    #     wavelengths = [696.54, 706.72]
    #     lines = [line1, line2]
    #     #2nd-order polynomial
    #     p2 = np.poly1d(np.polyfit(lines, wavelengths, 2))
    # else:
    #     print "No HgAr lamp calibration images found"
    # felist = glob.glob('*Fe*.fits')
    # if felist != []:
    #     fe = medcomb(felist)
    # else:
    #     print "No Fe lamp calibration images found"
