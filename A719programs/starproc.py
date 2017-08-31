#cmd line: 'python addframes.py filename top%'

import pyfits
import sys
import numpy as np
from scipy.stats import threshold
import matplotlib.pyplot as plt
import heapq

objname = sys.argv[1].split('.')[0]

#This loop is for the SNR vs % of images plot
snr = []
strl = []
for i in np.linspace(1, 100, 100):

	img = pyfits.open(sys.argv[1])

	#pct = float(sys.argv[2])
	pct = i

	#read pixel scale from header
	pxscl = float(img[0].header['PIXELSC'])

	#set diffraction limit ratio
	dlr = 0.015

	#find max point of first frame as reference
	firstmax = np.where(img[0].data==img[0].data.max())

	#define Strehl ratio equation
	def strehl(img):
		return (img.max() / np.sum(np.array(img))) / dlr

	#create array to store processed frames
	images = []
	#create array to store Strehl ratio values for each image
	strehl_frames = []
	for i in range(len(img)):
		img_data = img[i].data

		#choose region for background, subtract mean value
		bg_reg = img_data[:10, :10]
		bg_avg = np.mean(bg_reg, dtype=np.float64)
		img_data -= bg_avg

		#find max point of frame and shift to reference
		imax = np.where(img_data==img_data.max())
		xshift = int(np.average(firstmax[1])) - int(np.average(imax[1]))
		yshift = int(np.average(firstmax[0])) - int(np.average(imax[0]))
		img_data = np.roll(img_data, xshift, axis=1)
		img_data = np.roll(img_data, yshift, axis=0)	

		#find Strehl ratio of frame
		strehl_frames.append(strehl(img_data))

		#store processed frame
		images.append(img_data)

	#find top % of images
	top = heapq.nlargest(int(len(strehl_frames)*pct/100.0), strehl_frames)
	topindices = [i for i, j in enumerate(strehl_frames) if j in top]

	#select top % of images
	selectedimages = []
	for i in topindices:
		selectedimages.append(images[i])

	#save averaged image
	image_array = np.array(selectedimages)
#	if pct != 100:
#		append = '-lucky'+str(int(pct))
#	else:
#		append = ''
	img_avg = np.average(image_array, axis = 0)
#	hdu_avg = pyfits.PrimaryHDU(img_avg)
#	hdu_avg.writeto(objname+'_avg'+append+'.fits', clobber=True)

#calculate FWHM
	halfmax = img_avg.max()/2.0 #find value of half maximum
	red = threshold(img_avg, halfmax) #remove all pixels below half max
	profilex = red[firstmax[0], :] != 0 #find how many pixels above half max remain in row containing max pixel
	profiley = red[:, firstmax[1]] != 0 #find how many pixels above half max remain in column containing max pixel
	FWHM = pxscl * np.average([np.count_nonzero(profilex == True), np.count_nonzero(profiley == True)])
#	print 'FWHM =', FWHM

#calculate Strehl ratio of averaged image
	sr = strehl(img_avg)
	strl.append(sr)
#	print 'Strehl ratio =', sr

#plot Strehl ratios of individual frames as histogram
#	plt.figure()
#	plt.hist(strehl_frames, bins = len(strehl_frames)/4)
#	plt.title('Strehl ratios of individual frames for '+objname)
#	plt.xlabel('Strehl ratio')
#	plt.ylabel('# of frames')
#	plt.show()

	y, x = np.ogrid[:img[0].data.shape[0],:img[0].data.shape[1]]
	sigindex = (y-firstmax[0])**2.0 + (x-firstmax[1])**2 <= 4.0
	bgindex = (y-3)**2.0 + (x-3)**2.0 <= 4.0
	image_snr = np.sum(img_avg[sigindex])/(np.std(img_avg[bgindex]))

	snr.append(image_snr)

Calculate autocorrelation
unbiased = strehl_frames-np.mean(strehl_frames)
norm = np.sum(unbiased**2)
acor = np.correlate(unbiased, unbiased, "same")/norm

plt.figure()
plt.plot(acor)
plt.title('Autocorrelation of Strehl ratios by frame for '+objname)
plt.xlabel('Frame #')
plt.ylabel('Autocorrelation')
plt.show()

plt.figure()
plt.plot(np.linspace(1, 100, 100), snr)
plt.title('SNR for compound image using N% of frames for '+objname)
plt.xlabel('% of frames')
plt.ylabel('SNR')
plt.show()

plt.figure()
plt.plot(strl, snr)
plt.title('SNR vs. Strehl Ratio increasing % of images used for '+objname)
plt.xlabel('Strehl ratio with increasing % of images used')
plt.ylabel('SNR')
plt.show()
