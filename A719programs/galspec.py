#cmd line: python imcombine.py

import numpy as np
import sys
import pyfits
import glob
import matplotlib.pyplot as plt
from collections import Counter
from scipy.stats import threshold

#Median combine biases, calibrate and median combine images
biasnames = glob.glob('*.ZERO.fits')
biases = []
for bias in biasnames:
	bs = pyfits.open(bias)
	biases.append(bs[0].data)
bias_array = np.array(biases)
bias_med = np.median(bias_array, axis = 0)
imnames = glob.glob('*.rs0625bspec.fits')
images = []
im1 = pyfits.open(imnames[0])
pxscl = 0.15 #"/pix - googled
for img in imnames:
	im = pyfits.open(img)
	imdata = im[0].data - bias_med
	images.append(imdata)
image_array = np.array(images)
img_med = np.median(image_array, axis = 0)
img_med = threshold(img_med, threshmin = 0., threshmax = 700.)
hdu_med = pyfits.PrimaryHDU(img_med)
hdu_med.writeto('rs0625bspec_med.fits', clobber=True)

#Find midpoint of galaxy
mid = []
for i in range(np.shape(img_med)[1]):
	m = np.where(img_med[:, i] == img_med[1:np.shape(img_med)[1]-1, i].max())[0][0] #Necessary to cut off one pixel on each edge because edge pixels have strange values
	mid.append(m)
ymid = Counter(mid).most_common(1)[0][0]

#Find max Y across all X
x = np.arange(np.shape(img_med)[1])
y = []
for i in range(np.shape(img_med)[1]):
	ymax = np.where(img_med[:, i] == img_med[ymid-5:ymid+5, i].max())[0][0]
	y.append(ymax)

#3rd-order polynomial
p3 = np.poly1d(np.polyfit(x, y, 3))

#plt.figure()
#plt.ylim(ymid-5, ymid+5);
#plt.xlim(0, len(x));
#plt.xlabel('X pixel');
#plt.ylabel('Y pixel');
#plt.title('Location of max Y pixel across all X');
#plt.plot(x, y, '*', label = 'Max Pixel');
#plt.plot(p3(x), label = 'Polynomial Fit');
#plt.plot([0,len(x)], [ymid, ymid], '-', linewidth = 2.0, label = 'Galaxy Center');
#plt.legend();
#plt.show();

#Extract spectrum
spectrum = []
sub = []
for i in range(np.shape(img_med)[1]):
	lt = np.sum(img_med[int(p3(i)-9):int(p3(i)+9), i])
	spectrum.append(lt)
	sk = np.sum(img_med[int(p3(i)+373):int(p3(i)+391), i])
	sub.append(lt - sk)

#plt.figure()
#plt.xlabel('X Pixel');
#plt.ylabel('Intensity');
#plt.xlim(0, len(x));
#plt.title('Spectrum of Galaxy (raw)');
#plt.plot(x, spectrum);
#plt.show();

#plt.figure()
#plt.xlabel('X Pixel');
#plt.ylabel('Intensity');
#plt.xlim(0, len(x));
#plt.title('Spectrum of Galaxy (sky subtracted)');
#plt.plot(x, sub);
#plt.show();

#Extract lamp spectrum
lamp = pyfits.open('0317.HgAr.fits')[0].data
lampspectrum = []
for i in range(np.shape(lamp)[1]):
	lmp = np.sum(lamp[int(p3(i)-9):int(p3(i)+9), i])
	lampspectrum.append(lmp)

#Find line locations
lx = np.arange(np.shape(lamp)[1])
max1 = np.where(np.array(lampspectrum) == np.array(lampspectrum[:200]).max())[0][0]
line1 = np.sum(lampspectrum[max1-4:max1+4]*lx[max1-4:max1+4])/np.sum(lampspectrum[max1-4:max1+4])
max2 = np.where(np.array(lampspectrum) == np.array(lampspectrum[600:800]).max())[0][0]
line2 = np.sum(lampspectrum[max2-4:max2+4]*lx[max2-4:max2+4])/np.sum(lampspectrum[max2-4:max2+4])
max3 = np.where(np.array(lampspectrum) == np.array(lampspectrum[800:890]).max())[0][0]
line3 = np.sum(lampspectrum[max3-4:max3+4]*lx[max3-4:max3+4])/np.sum(lampspectrum[max3-4:max3+4])
max4 = np.where(np.array(lampspectrum) == np.array(lampspectrum[890:1000]).max())[0][0]
line4 = np.sum(lampspectrum[max4-4:max4+4]*lx[max4-4:max4+4])/np.sum(lampspectrum[max4-4:max4+4])
max5 = np.where(np.array(lampspectrum) == np.array(lampspectrum[1400:]).max())[0][0]
line5 = np.sum(lampspectrum[max5-4:max5+4]*lx[max5-4:max5+4])/np.sum(lampspectrum[max5-4:max5+4])
#print line1, line2, line3, line4, line5

#plt.figure()
#plt.xlabel('X Pixel');
#plt.ylabel('Intensity');
#plt.yscale('log');
#plt.annotate('435.83nm', xy = (max1-35, lampspectrum[max1]+30000));
#plt.annotate('546.07nm', xy = (max2-35, lampspectrum[max2]+40000));
#plt.annotate('576.90nm', xy = (max3-65, lampspectrum[max3]+5000));
#plt.annotate('579.06nm', xy = (max4, lampspectrum[max4]+7000));
#plt.annotate('696.54nm', xy = (max5-65, lampspectrum[max5]+50000));
#plt.xlim(0, len(lx));
#plt.title('Spectrum of HgAr Lamp');
#plt.plot(lx, lampspectrum);
#plt.show();

wavelengths = [435.83, 546.07, 576.90, 579.06, 696.54]
lines = [line1, line2, line3, line4, line5]

#2nd-order polynomial
p2 = np.poly1d(np.polyfit(lines, wavelengths, 2))

#Find Balmer line locations
bmax1 = np.where(np.array(sub) == np.array(sub[400:500]).max())[0][0]
balmer1 = np.sum(sub[bmax1-5:bmax1+5]*p2(x[bmax1-5:bmax1+5]))/np.sum(sub[bmax1-5:bmax1+5])
bmax2 = np.where(np.array(sub) == np.array(sub[1200:]).max())[0][0]
balmer2 = np.sum(sub[bmax2-10:bmax2+10]*p2(x[bmax2-10:bmax2+10]))/np.sum(sub[bmax2-10:bmax2+10])
bmax3 = np.where(np.array(sub) == np.array(sub[:250]).max())[0][0]
balmer3 = np.sum(sub[bmax3-5:bmax3+5]*p2(x[bmax3-5:bmax3+5]))/np.sum(sub[bmax3-5:bmax3+5])
#print balmer1, balmer2, balmer3

def z(obs, em):
	return (obs - em)/em

z1 = z(balmer1, 486.13)
z2 = z(balmer2, 656.28)
z3 = z(balmer3, 434.05)
#print z1, z2, z3

#plt.figure()
#plt.xlabel('Wavelength (nm)');
#plt.annotate(r'$H\beta$', xy = (p2(bmax1-10), sub[bmax1]+100));
#plt.annotate(r'$H\alpha$', xy = (p2(bmax2-10), sub[bmax2]+100));
#plt.annotate(r'$H\gamma$', xy = (p2(bmax3-10), sub[bmax3]+100));
#plt.ylabel('Intensity');
#plt.title('Spectrum of Galaxy');
#plt.plot(p2(x), sub);
#plt.show();

##Monte Carlo approx
def mcp2(a, b, c, x):
	return a*x**2 + b*x + c

N = 1000
mczHb = []
mczHa = []
mczHy = []
poly = np.zeros((3, N))	
Hapix = []
for n in range(N):
	mclpix1 = int(np.random.normal(max1, 4))
	mclamp1 = np.sum(lampspectrum[mclpix1-4:mclpix1+4]*lx[mclpix1-4:mclpix1+4])/np.sum(lampspectrum[mclpix1-4:mclpix1+4])
	mclpix2 = int(np.random.normal(max2, 4))
	mclamp2 = np.sum(lampspectrum[mclpix2-4:mclpix2+4]*lx[mclpix2-4:mclpix2+4])/np.sum(lampspectrum[mclpix2-4:mclpix2+4])
	mclpix3 = int(np.random.normal(max3, 4))
	mclamp3 = np.sum(lampspectrum[mclpix3-4:mclpix3+4]*lx[mclpix3-4:mclpix3+4])/np.sum(lampspectrum[mclpix3-4:mclpix3+4])
	mclpix4 = int(np.random.normal(max4, 4))
	mclamp4 = np.sum(lampspectrum[mclpix4-4:mclpix4+4]*lx[mclpix4-4:mclpix4+4])/np.sum(lampspectrum[mclpix4-4:mclpix4+4])
	mclpix5 = int(np.random.normal(max5, 4))
	mclamp5 = np.sum(lampspectrum[mclpix5-4:mclpix5+4]*lx[mclpix5-4:mclpix5+4])/np.sum(lampspectrum[mclpix5-4:mclpix5+4])
	mclines = [mclamp1, mclamp2, mclamp3, mclamp4, mclamp5]
	mcp = np.polyfit(mclines, wavelengths, 2)
	poly[0, n] = mcp[0]
	poly[1, n] = mcp[1]
	poly[2, n] = mcp[2]
	mcpix1 = int(np.random.normal(bmax1, 5))
	mcbalmer1 = np.sum(sub[mcpix1-5:mcpix1+5]*mcp2(poly[0, n], poly[1, n], poly[2, n], x[mcpix1-5:mcpix1+5]))/np.sum(sub[mcpix1-5:mcpix1+5])
	mcz1 = z(mcbalmer1, 486.13)
	mczHb.append(mcz1)
	mcpix2 = int(np.random.normal(bmax2, 10))
	Hapix.append(mcpix2)
	mcbalmer2 = np.sum(sub[mcpix2-10:mcpix2+10]*mcp2(poly[0, n], poly[1, n], poly[2, n], x[mcpix2-10:mcpix2+10]))/np.sum(sub[mcpix2-10:mcpix2+10])
	mcz2 = z(mcbalmer2, 656.28)	
	mczHa.append(mcz2)
	mcpix3 = int(np.random.normal(bmax3, 5))
	mcbalmer3 = np.sum(sub[mcpix3-5:mcpix3+5]*mcp2(poly[0, n], poly[1, n], poly[2, n], x[mcpix3-5:mcpix3+5]))/np.sum(sub[mcpix3-5:mcpix3+5])
	mcz3 = z(mcbalmer3, 434.05)
	mczHy.append(mcz3)

zHb = float('%s' % float('%.3g' % np.average(mczHb)))
sigzHb = float('%s' % float('%.2g' % np.std(mczHb)))
zHa = float('%s' % float('%.3g' % np.average(mczHa)))
sigzHa = float('%s' % float('%.2g' % np.std(mczHa)))
zHy = float('%s' % float('%.3g' % np.average(mczHy)))
sigzHy = float('%s' % float('%.2g' % np.std(mczHy)))

mcz = mczHb + mczHa + mczHy
totz = float('%s' % float('%.3g' % np.average(mcz)))
sigz = float('%s' % float('%.2g' % np.std(mcz)))

plt.figure()
plt.title('Histogram of possible redshift values for each Balmer line');
plt.xlabel('Redshift (z)');
plt.ylabel('Instances');
plt.hist(mczHa, 100, label=r'$H\alpha; z = %s \pm %s$'%(zHa, sigzHa));
plt.hist(mczHb, 100, label=r'$H\beta; z = %s \pm %s$'%(zHb, sigzHb));
plt.hist(mczHy, 100, label=r'$H\gamma; z = %s \pm %s$'%(zHy, sigzHy));
plt.legend();
plt.show();

plt.figure()
plt.title('Histogram of possible redshift values for combined data; $z = %s \pm %s$'%(totz, sigz));
plt.xlabel('Redshift (z)');
plt.ylabel('Instances');
plt.hist(mcz, 100);
plt.show();

v = []
verr = []
loc = []
wave = []
cent = -100
while cent < 100:
	loc.append(cent*pxscl)
	mcv = []
	spec = []
	w = []
	for i in range(np.shape(img_med)[1]):
		lt = np.sum(img_med[int((p3(i)+cent)-3):int((p3(i)+cent)+3), i])
		sk = np.sum(img_med[int(p3(i)+373):int(p3(i)+391), i])
		spec.append(lt-sk)
	for n in range(N):
		mcbalmer2 = np.sum(spec[Hapix[n]-10:Hapix[n]+10]*mcp2(poly[0, n], poly[1, n], poly[2, n], x[Hapix[n]-10:Hapix[n]+10]))/np.sum(spec[Hapix[n]-10:Hapix[n]+10])
		w.append(mcbalmer2)
		mv = (z(mcbalmer2, 656.28) - mczHa[n])*299792.458 #c in km/s
		mcv.append(mv)
	wave.append(np.average(w))	
	v.append(np.average(mcv))
	verr.append(np.std(mcv))
	cent += 1

plt.figure()
plt.xlabel('Angluar shift from galaxy center (arcseconds)')
plt.ylabel('Rotational velocity (km/s)')
plt.title('Rotation curve of galaxy')
plt.errorbar(loc, v, yerr=verr);
plt.show();

plt.figure()
plt.xlabel('Angluar shift from galaxy center (arcseconds)')
plt.ylabel('Rotational velocity (km/s)')
plt.title('Rotation curve of galaxy (no error bars to see detail)')
plt.plot(loc, v);
plt.show();

plt.figure()
plt.xlabel('Angluar shift from galaxy center (arcseconds)')
plt.ylabel(r'Wavelength of $H_\alpha (nm)$')
plt.title('Wavelength vs Offset from Galaxy')
plt.plot(loc, wave);
plt.show();
