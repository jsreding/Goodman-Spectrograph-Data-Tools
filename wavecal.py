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

def gauss(r, p):
    p = int(p)
    r = int(r)
    l = GaussianModel(prefix='l_')
    pars = l.guess(imdata[r, p-5:p+5], x=x[p-5:p+5])
    pars.update(l.make_params())
    pars['l_center'].set(p, min=p-5, max=p+5)
    pars['l_sigma'].set(0.5, max=np.sqrt(3))
    pars['l_amplitude'].set(50, min=0)
    mod = l
    # plt.plot(x[p-5:p+5], imdata[r, p-5:p+5], label='data')
    out = mod.fit(imdata[r, p-5:p+5], pars, x=x[p-5:p+5])
    comps = out.eval_components(x=x[p-5:p+5])
    # print(out.fit_report(min_correl=0.5))
    # plt.plot(x[p-5:p+5], out.best_fit, 'r-', label='best fit')
    # plt.plot(x[p-5:p+5], comps['l_'], 'b--', label='line')
    # plt.legend()
    # plt.show()
    return out.params['l_center']
###
imdata = pyfits.open('0018_Fe_930_blue_long.fits')[0].data

#Note: Cut off first and last 50 pixels for readout data
biasnames = glob.glob('*Zero*.fits')
bias = medcomb(biasnames)
flatnames = glob.glob('*QuFlat*.fits')
flat = medcomb(flatnames)
if np.shape(imdata)[1] != np.shape(bias)[1]:
    bias = ndimage.zoom(bias, (1, np.shape(imdata)[1]/np.shape(bias)[1]), order=0)[:,50:-50]
if np.shape(imdata)[1] != np.shape(flat)[1]:
    flat = ndimage.zoom(flat, (1, np.shape(imdata)[1]/np.shape(flat)[1]), order=0)[:,50:-50]

# To possibly change - use only arcs taken right before or after science image

imdata = (imdata[:,50:-50] - bias)/(flat - bias) * np.average(flat-bias)
# hdu = pyfits.PrimaryHDU(imdata)
# hdu.writeto('clean_0018.fits', clobber=True)

# Choose your setup
setup = raw_input("Setup [930B, 930R, 1200Ha, 1200Ca]: ")
if setup == '930B':
    gr, grang, camang = 930., 12., 24.
elif setup == '930R':
    gr, grang, camang = 930., 20., 35.2
elif setup == '1200Ha':
    gr, grang, camang = 1200., 24., 48.
elif setup == '1200Ca':
    gr, grang, camang = 1200., 31., 62.
#Central wavelength
print "Grating angle:", grang
print "Camera angle:", camang
cwv = np.sin(grang*np.pi/180)*20000000/gr
print "Central wavelength (A):", np.round(cwv, 3)

#From Josh Fuchs; ([930B], [930R])
#Fe line typical pixel locations (adjusted for readout cutoff)
fepxsb = np.array([431.795, 451.264, 942.966, 1057.76, 1174.6, 1194.97, 1315.35, 1381.1, 1444.58, 1457.81, 1538.65, 1544.11, 1630.61, 1682.99, 1713.34, 1726.03, 1779.61, 1886.38, 1893.19, 1959.28, 1968.19, 1980.81, 2018.27, 2078.23, 2088.47, 2132.53, 2194.03, 2210.85, 2279.64, 2361.34, 2443.06, 2468.22, 2515.14, 2630.53, 2795.45, 2807.86, 2817.11, 2886.45, 2985.15, 3085.52, 3162.56, 3184.41, 3367.86, 3493.34, 3602.21, 3795.76, 3845.65, 3907.57])-50.
fepxsr = np.array([1155.1, 102.964, 142.88, 276.362, 438.35, 532.819, 631.475, 755.5, 798.185, 831.719, 1062.43, 1086.89, 1249.02, 1316.94, 1475.64, 1566.37, 1762.42, 1910, 2053.55, 2072.75, 2168.64, 2179.14, 2271.52, 2318.23, 2347.65, 2370.04, 2417.91, 2645.82, 2672.71, 3385.8, 3562.11, 3620.55, 3765.62, 3913.77, 3935.87, 4016.2, 4034.6])-50.

#Fe line wavelength locations
fewvsb = np.array([3729.3087, 3737.1313, 3946.0971, 3994.7918, 4044.4179, 4052.9208, 4103.9121, 4131.7235, 4158.5905, 4164.1795, 4198.3036, 4200.6745, 4237.2198, 4259.3619, 4271.7593, 4277.5282, 4300.1008, 4345.168, 4348.064, 4375.9294, 4379.6668, 4385.0566, 4400.9863, 4426.0011, 4430.189, 4448.8792, 4474.7594, 4481.8107, 4510.7332, 4545.0519, 4579.3495, 4589.8978, 4609.5673, 4657.9012, 4726.8683, 4732.0532, 4735.9058, 4764.8646, 4806.0205, 4847.8095, 4879.8635, 4889.0422, 4965.0795, 5017.1628, 5062.0371, 5141.7827, 5162.2846, 5187.7462])
fewvsr = np.array([6043.223, 6466.5526, 6483.0825, 6538.112, 6604.8534, 6643.6976, 5875.618, 6684.2929, 6752.8335, 6766.6117, 6861.2688, 6871.2891, 6937.6642, 6965.4307, 7030.2514, 7067.2181, 7147.0416, 7206.9804, 7265.1724, 7272.9359, 7311.7159, 7316.005, 7353.293, 7372.1184, 7383.9805, 7392.9801, 7412.3368, 7503.8691, 7514.6518, 7798.5604, 7868.1946, 7891.075, 7948.1964, 8006.1567, 8014.7857, 8046.1169, 8053.3085])

x = np.arange(np.shape(imdata)[1])
if setup == '930B':
    pxs, wvs = fepxsb, fewvsb
elif setup == '930R':
    pxs, wvs = fepxsr, fewvsr
lnslopes = np.zeros(len(pxs))
for l in range(len(pxs)):
    n = 0
    lc = np.zeros(20)
    for y in np.linspace(10, 199, 20):
        lc[n] = gauss(y, pxs[l])
        n += 1
    lnslopes[l] = np.polyfit(np.linspace(10, 199, 20), lc, 1)[0]
print np.mean(lnslopes)
