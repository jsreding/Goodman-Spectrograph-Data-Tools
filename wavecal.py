import numpy as np
from scipy import ndimage
import os
import pyfits
import glob
from astropy import stats
import matplotlib.pyplot as plt
import lacosmic
from collections import Counter
from lmfit.models import GaussianModel

###
def medcomb(names):
    ims = []
    for n in range(len(names)):
        ims.append(pyfits.open(names[n])[0].data)
    im_med = np.median(np.array(ims), axis=0)
    return im_med

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
###
# Choose your setup
setup = raw_input("Setup [930B, 930R, 1200Ha, 1200Ca]: ")
if setup == '930B':
    gr, grang, camang = 930., 12., 24.
    arcnames = glob.glob("*blue_long.fits")
elif setup == '930R':
    gr, grang, camang = 930., 20., 35.2
    arcnames = glob.glob("*red_long.fits")
elif setup == '1200Ha':
    gr, grang, camang = 1200., 24., 48.
elif setup == '1200Ca':
    gr, grang, camang = 1200., 31., 62.
else:
    print "Setup not recognized"
    quit()
#Central wavelength
print "Grating angle:", grang
print "Camera angle:", camang
cwv = np.sin(grang*np.pi/180)*20000000/gr
print "Central wavelength (A):", np.round(cwv, 3)

arcdata = medcomb(arcnames)

#Note: Cut off first and last 50 pixels for readout data
biasnames = glob.glob('*Zero*.fits')
bias = medcomb(biasnames)
flatnames = glob.glob('*Q*Flat*.fits')
flat = medcomb(flatnames)
if np.shape(arcdata)[1] != np.shape(bias)[1]:
    bias = ndimage.zoom(bias, (1, np.shape(arcdata)[1]/np.shape(bias)[1]), order=0)[:,50:-50]
if np.shape(arcdata)[1] != np.shape(flat)[1]:
    flat = ndimage.zoom(flat, (1, np.shape(arcdata)[1]/np.shape(flat)[1]), order=0)[:,50:-50]

# To possibly change - use only arcs taken right before or after science image

arcdata = (arcdata[:,50:-50] - bias)/(flat - bias) * np.average(flat-bias)
# hdu = pyfits.PrimaryHDU(arcdata)
# hdu.writeto('clean_blue.fits', clobber=True)

#Fe line typical pixel locations and wavelengths, by setup (adjusted for readout cutoff)
shift = np.argmax(arcdata[100,1257:1357])-50
#This line is at least 50 pix from the next bright one
fepxsb = np.array([372., 379., 1121., 1307., 1370., 1384., 1464., 1470., 1609., 1638., 1652., 1705., 1784., 1818., 2204., 2285., 2439., 2554., 2718., 2739., 2808., 3084.])+shift
fewvsb = np.array([3734.8636, 3737.1313, 4045.8130, 4131.7235, 4158.5905, 4164.1795, 4198.3036, 4200.6745, 4259.3619, 4271.7593, 4277.5282, 4300.1008, 4333.5612, 4348.064, 4510.7332, 4545.0519, 4609.5673, 4657.9012, 4726.8683, 4735.9058, 4764.8648, 4879.8635])

x = np.arange(np.shape(arcdata)[1])
if setup == '930B':
    pxs, wvs = fepxsb, fewvsb
elif setup == '930R':
    pxs, wvs = fepxsr, fewvsr
lines = np.zeros((2, len(pxs)))
for l in range(len(pxs)):
    n = 0
    lc = np.zeros(20)
    for y in np.linspace(10, 199, 20):
        lc[n] = arcgauss(y, pxs[l])
        n += 1
    lines[0][l] = np.polyfit([int(p) for p in np.linspace(10, 199, 20)], lc, 1)[0]
    lines[1][l] = np.polyfit([int(p) for p in np.linspace(10, 199, 20)], lc, 1)[1]
#This value 'loc' will be the average y-position of the science object in the data image
for loc in np.linspace(80, 120, 5):
    pixlocs = lines[1] + loc*lines[0]
    wavesol = np.poly1d(np.polyfit(pixlocs, wvs, 3))+50
    print wavesol
