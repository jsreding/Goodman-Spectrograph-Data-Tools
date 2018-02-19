### Goodman Spectrograph image reduction and spectrum extraction (with cosmic ray removal)  ###
### Joshua Reding UNC-CH 7/27/17   ###

import numpy as np
import os
import pyfits
import glob
from astropy import stats
from scipy import ndimage
import matplotlib.pyplot as plt
import lacosmic
from collections import Counter
from lmfit.models import GaussianModel

###################################################
def object():
    objname = raw_input("Object name: ")
    setup = raw_input("Setup [930_blue, 930_red]: ")
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
    data = np.float64(datalist[0].data)
    gain = hdr['GAIN'] #1.33 from 2017-06-07
    rdnoise = hdr['RDNOISE']
    c = lacosmic.lacosmic(data, contrast=2, cr_threshold=int(crt), neighbor_threshold=.5, readnoise=rdnoise, effective_gain=gain, maxiter=4)
    return c[0][1:-1, 50:], hdr

def arcgauss(r, p):
    p = int(p)
    r = int(r)
    l = GaussianModel(prefix='l_')
    pars = l.guess(arcdata[r, p-5:p+5], x=xarc[p-5:p+5])
    pars.update(l.make_params())
    pars['l_center'].set(p, min=p-5, max=p+5)
    pars['l_sigma'].set(0.5, max=np.sqrt(3))
    pars['l_amplitude'].set(50, min=0)
    mod = l
    # plt.plot(xarc[p-5:p+5], arcdata[r, p-5:p+5], label='data')
    out = mod.fit(arcdata[r, p-5:p+5], pars, x=xarc[p-5:p+5])
    comps = out.eval_components(x=xarc[p-5:p+5])
    # print(out.fit_report(min_correl=0.5))
    # plt.plot(xarc[p-5:p+5], out.best_fit, 'r-', label='best fit')
    # plt.plot(xarc[p-5:p+5], comps['l_'], 'b--', label='line')
    # plt.legend()
    # plt.show()
    return out.params['l_center']
###################################################
oldfiles = glob.glob('clean*.fits')+glob.glob('spec*.fits')
for o in oldfiles:
    os.remove(o)

objname, setup = object()
print "Finding spectra for "+objname+" taken in "+setup+"setup..."

if setup == '930_blue':
    gr, grang, camang = 930., 12., 24.
    arcnames = glob.glob("*blue_long.fits")
#####Soon to come#####
# elif setup == '930_red':
#     gr, grang, camang = 930., 20., 35.2
#     arcnames = glob.glob("*red_long.fits")
# elif setup == '1200Ha':
#     gr, grang, camang = 1200., 24., 48.
# elif setup == '1200Ca':
#     gr, grang, camang = 1200., 31., 62.
else:
    print "Setup not recognized"
    quit()

imnames = glob.glob('*'+objname+'_'+setup+'.fits')
imnames = sorted(imnames, key=lambda imsa: int(imsa.split('_')[0]))

biasnames = glob.glob('*Zero*.fits')
bias_med = medcomb(biasnames)[1:-1, 50:]
print np.shape(bias_med)
flatnames = glob.glob('*Flat*'+setup+'.fits')
flat_med = medcomb(flatnames)[1:-1:, 50:]
print np.shape(flat_med)

print "Grating angle:", grang
print "Camera angle:", camang
cwv = np.sin(grang*np.pi/180)*20000000/gr
print "Central wavelength (A):", np.round(cwv, 3)

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
    	y[i] = ymid-7 + np.argmax(imdata[ymid-7:ymid+7, i])
    #3rd-order polynomial
    sp = np.poly1d(np.polyfit(x, y, 3))

    #Extract spectrum
    spectrum = np.zeros(np.shape(imdata)[1])
    for i in range(np.shape(imdata)[1]):
            x = np.array(np.ndarray.tolist(np.arange(20))+np.ndarray.tolist(np.arange(40,60)))
            y = np.array(np.ndarray.tolist(imdata[int(sp(i)-30):int(sp(i)-10), i]) + np.ndarray.tolist(imdata[int(sp(i)+10):int(sp(i)+30), i]))
            sk = np.poly1d(np.polyfit(x, y, 3))
            spectrum[i] = np.sum(imdata[int(sp(i)-10):int(sp(i)+10), i])
            # spectrum[i] = np.sum(imdata[int(sp(i)-10):int(sp(i)+10), i] - sk(np.arange(20,40)))
    spectrum = stats.sigma_clip(spectrum, sigma=5)

    ###Wavelength Calibration###
    arcdata = medcomb(arcnames)[1:-1, 50:-50]

    if np.shape(arcdata)[1] != np.shape(bias_med)[1]:
        bias_med = ndimage.zoom(bias_med, (1, np.shape(arcdata)[1]/np.shape(bias_med)[1]), order=0)
    if np.shape(arcdata)[1] != np.shape(flat_med)[1]:
        flat_med = ndimage.zoom(flat_med, (1, np.shape(arcdata)[1]/np.shape(flat_med)[1]), order=0)

    # To possibly change - use only arcs taken right before or after science image

    arcdata = (arcdata - bias_med)/(flat_med - bias_med) * np.average(flat_med-bias_med)
    # hdu = pyfits.PrimaryHDU(arcdata)
    # hdu.writeto('clean_blue.fits', clobber=True)

    xarc = np.arange(np.shape(arcdata)[1])
    #This line is at least 50 pix from the next bright one
    shift = int(arcgauss(100, 1257+np.argmax(arcdata[100,1257:1357]))-1307)
    #Fe line typical pixel locations and wavelengths, by setup (adjusted for readout cutoff)
    fepxsb = np.array([1121., 1307., 1370., 1384., 1609., 1638., 1652., 1705., 1784., 1818., 2204., 2285., 2439., 2739., 2808., 3084.])+shift
    fewvsb = np.array([4045.8130, 4131.7235, 4158.5905, 4164.1795, 4259.3619, 4271.7593, 4277.5282, 4300.1008, 4333.5612, 4348.064, 4510.7332, 4545.0519, 4609.5673, 4735.9058, 4764.8648, 4879.8635])

    if setup == '930_blue':
        pxs, wvs = fepxsb, fewvsb
    # elif setup == '930_red':
    #     pxs, wvs = fepxsr, fewvsr
    lines = np.zeros((2, len(pxs)))
    for l in range(len(pxs)):
        n = 0
        lc = np.zeros(20)
        for v in np.linspace(10, 190, 20):
            lc[n] = arcgauss(v, pxs[l])
            n += 1
        lines[0][l] = np.polyfit([int(p) for p in np.linspace(10, 199, 20)], lc, 1)[0]
        lines[1][l] = np.polyfit([int(p) for p in np.linspace(10, 199, 20)], lc, 1)[1]
    pixlocs = lines[1] + np.argmax(imdata[:, 100])*lines[0]
    wavesol = np.poly1d(np.polyfit(pixlocs, wvs, 3))
    ###

    plt.figure()
    plt.title(img)
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Flux (counts)')
    plt.plot(wavesol(np.arange(len(xscale))), spectrum)
    plt.show()

# plt.figure()
# plt.scatter(mcl, wvs, label='lines')
# plt.plot(x, wvx, label='fit')
# plt.legend()
# plt.show()
# print wvs - mcpoly(mcl)
#
# hdu = pyfits.PrimaryHDU(final)
# hdu.header = hdr
# hdu.header['BZERO'] = 0
# hdu.writeto('spec_'+img.split('.')[0]+'.fits', clobber=True)
