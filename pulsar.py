import numpy as np
from scipy.ndimage.filters import gaussian_filter
import pyfits
import glob
import matplotlib.pyplot as plt

obs = glob.glob('acisf10135N003_evt2.fits') # Just one image
#obs = glob.glob('acisf01994N004_evt2.fits') # Crab pulsar
#obs = glob.glob('acis*.fits') # All images for movie
#obs_sorted = sorted(obs, key=lambda obs: obs.split('N')[0])
grades = [0, 2, 4, 6]
for o in obs:
#for o in obs_sorted:
	im = pyfits.open(o)
	im_status = np.array(im[1].data['STATUS'])
	im_grade = np.array(im[1].data['GRADE'])
	im_energy = np.array(im[1].data['ENERGY'])
	events = np.arange(len(im_status))
	failed = []
	for e in events:
		if np.sum(im_status[e, :]) != 0 or (im_grade[e] in grades) == False or im_energy[e] < 1200:
			failed.append(e)
	events = np.delete(events, failed)
	im_X = np.divide(np.array(im[1].data['X'][events]), 4.).astype(int)
	im_Y = np.divide(np.array(im[1].data['Y'][events]), 4.).astype(int)
	im_energy = np.array(im[1].data['ENERGY'][events])
	frame = np.zeros((im_Y.max(), im_X.max()))
	eframe = np.zeros((im_Y.max(), im_X.max()))
	for i in range(len(im_X)):
		frame[im_Y[i]-1, im_X[i]-1] += 1
		eframe[im_Y[i]-1, im_X[i]-1] += (im_energy[i]/1000.)
	smooth = gaussian_filter(frame, sigma = 2)
	hdu = pyfits.PrimaryHDU(frame)
	hdu.writeto('nobg-'+o, clobber=True)
	hdu = pyfits.PrimaryHDU(smooth)
	hdu.writeto('sm_nobg-'+o, clobber=True)

#	if o == 'acisf10135N003_evt2.fits':
#		calmax = np.where(smooth == smooth.max())
##		hdu = pyfits.PrimaryHDU(smooth)
##		hdu.writeto('mov-'+o, clobber=True)

#	jetmax = [1082, 1017] # Manually estimated
#	R = 10
#	y, x = np.ogrid[:smooth.shape[0],:smooth.shape[1]]
#	pulsaraperture = (y-calmax[0][0])**2.0 + (x-calmax[1][0])**2 <= R**2
#	jetaperture = (y-jetmax[0])**2.0 + (x-jetmax[1])**2 <= R**2
#	
#	# Because there were large outliers
#	pulsar = []	
#	for i in eframe[pulsaraperture]:
#		if i < 1500.:
#			pulsar.append(i)
#	
#	plt.figure()
#	plt.title('Histogram of energies for pulsar aperture');
#	plt.xlabel('Energy (KeV)');
#	plt.ylabel('Instances');
#	plt.hist(pulsar, 50);
#	plt.show();
#	
#	plt.figure()
#	plt.title('Histogram of energies for jet aperture');
#	plt.xlabel('Energy (KeV)');
#	plt.ylabel('Instances');
#	plt.hist(eframe[jetaperture], 50);
#	plt.show();

	# Align other images to make movie
#	else:
#		m = np.where(smooth == smooth.max())
#		xshift = calmax[1][0] - m[1][0]
#		yshift = calmax[0][0] - m[0][0]
#		smooth = np.roll(smooth, xshift, axis=1)
#		smooth = np.roll(smooth, yshift, axis=0)
#		hdu = pyfits.PrimaryHDU(smooth)
#		hdu.writeto('mov-'+o, clobber=True)

