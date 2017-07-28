### QUICK TRANSFER CURVE PRODUCTION TOOL###
#####	Joshua Reding - 06/09/17	#####

# Naming scheme must be *_transfcurve_*_freq_attn.fits for flats and *ZERO*freq_attn.fits for biases

import pyfits
import numpy as np
import matplotlib.pyplot as plt
import re
import glob

files = glob.glob('*_transfcurve_*_*_*.fits')
freq = set([x.split('_')[3] for x in files]) #Tabulate frequencies used from filenames
attn = set([re.split('[_ .]', x)[4] for x in files]) #Tabulate attentuations used from filenames
avgs = []
vrcs = []
f = '200'
a = '0'
#for f in freq:
#	for a in attn:
centers = [[250,700],[300,850],[550,725],[750,700],[850,850],[900,725],[1150,750],[1200,850],[1300,725]]
for choice in centers:
        xcenter = choice[0]
        ycenter = choice[1]
        imnames = glob.glob('*transfcurve*'+f+'_'+a+'.fits')
        imnames_sorted = sorted(imnames, key=lambda imnames: int(imnames.split('_')[0]))
        biasnames = glob.glob('*ZERO*'+f+'_'+a+'.fits')
        biases = []
        avgs = []
        vrcs = []
        #Median combine biases
        for b in biasnames:
                bias = pyfits.open(b)
                biases.append(bias[0].data)
        biasarray = np.array(biases)
        biasmed = np.median(biasarray)
        #Process exposure pairs
        n = 0
        while n <= len(imnames)/2:
                img1 = pyfits.open(imnames_sorted[n])
                img2 = pyfits.open(imnames_sorted[n+1])
                if img1[0].header['PG4_0'] == img2[0].header['PG4_0']:
                        im1data = img1[0].data
                        im2data = img2[0].data
                        diff = im1data - im2data
                        #ycenter = int(np.shape(diff)[0]/2.)
                        #xcenter = int(np.shape(diff)[1]/2.)
                        #print xcenter, ycenter
                        avg = np.mean(im1data[ycenter-50:ycenter+50, xcenter-50:xcenter+50])-biasmed
                        vrc = np.var(diff[ycenter-50:ycenter+50, xcenter-50:xcenter+50])/2.
                        avgs.append(avg)
                        vrcs.append(vrc)
                        n += 2
                else:
                        n += 1
        #Best fit line (slope is reciprocal of gain)
        poly = np.polyfit(avgs[0:8], vrcs[0:8], 1)
        gain = 1./poly[0]
        p = np.poly1d(poly)
        print choice
        print gain
        print ''
        '''
        x = np.linspace(0, 60000, 1000)
        plt.figure()
        plt.title('Transfer curve for '+f+' kHz, ATTN '+a)
        plt.xlabel('Average counts in 100x100 box at center')
        plt.ylabel('Variance of 100x100 box at center')
        plt.plot(x, p(x), label = 'Gain = %s'%(gain))
        plt.scatter(avgs, vrcs)
        plt.legend()
        plt.show()
        '''
