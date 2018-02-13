import numpy as np
import pyfits
import glob
import matplotlib.pyplot as plt

#read in bias
bsa = pyfits.open('Prompt5_mbias_351001.fits')
bsb = pyfits.open('Prompt5_mbias_350891.fits')

#read in dark
dk = pyfits.open('Prompt5_mdark_80.0_350907.fits')

#read in flat and reduce
ft = pyfits.open('Prompt5_mflat_None_350976.fits')

#Start loop to read filters
imsa = glob.glob('wasp_43_1967313_V_???.fits')
imsa_sort = sorted(imsa, key=lambda imsa: int((imsa.split('_')[4]).split('.')[0]))
imsb = glob.glob('wasp_43_cont_1968840_V_???.fits')
imsb_sort = sorted(imsb, key=lambda imsb: int((imsb.split('_')[5]).split('.')[0]))
sets = [imsa_sort, imsb_sort]
wasp = []
waspa = []
ref1 = []
ref2 = []
ref3 = []
ref4 = []
ref5 = []
ref1a = []
ref2a = []
ref3a = []
ref4a = []
ref5a = []
time = []
timea = []
for s in sets:
	for im in s:
		lt = pyfits.open(im)

		#Bias subtract dark, flat, and image by proper image set bias
		if s == imsa_sort:
			if ('BIASCORR' in dk[0].header) == False:
				dk_bs = (dk[0].data - bsa[0].data)
			else:
				dk_bs = dk[0].data
			if ('BIASCORR' in ft[0].header) == False:
				ft_bs = ft[0].data - bsa[0].data
			else:
				ft_bs = ft[0].data
			if ('BIASCORR' in lt[0].header) == False:
				lt_bs = lt[0].data - bsa[0].data
			else:
				lt_bs = lt[0].data
		else:
			if ('BIASCORR' in dk[0].header) == False:
				dk_bs = (dk[0].data - bsb[0].data)
			else:
				dk_bs = dk[0].data
			if ('BIASCORR' in ft[0].header) == False:
				ft_bs = ft[0].data - bsb[0].data
			else:
				ft_bs = ft[0].data
			if ('BIASCORR' in lt[0].header) == False:
				lt_bs = lt[0].data - bsb[0].data
			else:
				lt_bs = lt[0].data

		#Scale dark
		dk_red = dk_bs * (lt[0].header['EXPOSURE']/dk[0].header['EXPOSURE'])

		#Dark correct flat
		if ('DARKCORR' in ft[0].header) == False:
			ft_red = (ft_bs - dk_red)/np.average((ft_bs - dk_red))
		else:
			ft_red = ft_bs/np.average(ft_bs)

		#Dark correct image
		if ('DARKCORR' in lt[0].header) == False:
			lt_red = lt_bs - dk_red
		else:
			lt_red = lt_bs
		#Reduce image and save
		lt_data = (lt_red - dk_red) / ft_red
#		hdu = pyfits.PrimaryHDU(lt_data)
#		hdu.writeto('red-'+im, clobber=True)

		#Find WASP-43 and references
		waspy, waspx = np.where(lt_data == lt_data[495:595, 370:470].max())[0][0], np.where(lt_data == lt_data[495:595, 370:470].max())[1][0]
		ref1y, ref1x = np.where(lt_data == lt_data[610:710, 880:980].max())[0][0], np.where(lt_data == lt_data[610:710, 880:980].max())[1][0]
		ref2y, ref2x = np.where(lt_data == lt_data[139:239, 784:884].max())[0][0], np.where(lt_data == lt_data[139:239, 784:884].max())[1][0]
		ref3y, ref3x = np.where(lt_data == lt_data[761:861, 187:287].max())[0][0], np.where(lt_data == lt_data[761:861, 187:287].max())[1][0]
		ref4y, ref4x = np.where(lt_data == lt_data[378:478, 516:616].max())[0][0], np.where(lt_data == lt_data[378:478, 516:616].max())[1][0]
		ref5y, ref5x = np.where(lt_data == lt_data[96:196, 358:458].max())[0][0], np.where(lt_data == lt_data[96:196, 358:458].max())[1][0]

		#Photometer and background subtract WASP-43 and reference
		R1 = 20
		R2 = 40
		y, x = np.ogrid[:lt_data.shape[0],:lt_data.shape[1]]

		waspaperture = (y-waspy)**2.0 + (x-waspx)**2 <= R1**2
		waspraw = lt_data[waspaperture]
		waspannulus = ((y-waspy)**2.0 + (x-waspx)**2 <= R2**2) - waspaperture 
		#BG sigma clipping		
		waspbg = np.median(lt_data[waspannulus][np.where(abs(lt_data[waspannulus]-np.average(lt_data[waspannulus]) < np.std(lt_data[waspannulus])))])
#		waspbg = np.median(lt_data[waspannulus])
		waspsub = [i - waspbg for i in waspraw]
		waspfinal = np.sum(waspsub)
		
		def photometer(rx, ry):
			refaperture = (y-ry)**2.0 + (x-rx)**2 <= R1**2
			refraw = lt_data[refaperture]
			refannulus = ((y-ry)**2.0 + (x-rx)**2 <= R2**2) - refaperture 
			#BG sigma clipping
			refbg = np.median(lt_data[refannulus][np.where(abs(lt_data[refannulus]-np.average(lt_data[refannulus]) < np.std(lt_data[refannulus])))])
#			refbg = np.median(lt_data[refannulus])
			refsub = [i - refbg for i in refraw]
			return np.sum(refsub)
			
		ref1final = photometer(ref1x, ref1y)
		ref2final = photometer(ref2x, ref2y)
		ref3final = photometer(ref3x, ref3y)
		ref4final = photometer(ref4x, ref4y)
		ref5final = photometer(ref5x, ref5y)
		
		if s == imsa_sort:
			waspa.append(waspfinal)
			ref1a.append(ref1final)
			ref2a.append(ref2final)
			ref3a.append(ref3final)
			ref4a.append(ref4final)
			ref5a.append(ref5final)
			timea.append(lt[0].header['JD'])
		else:
			wasp.append(waspfinal)
			ref1.append(ref1final)
			ref2.append(ref2final)
			ref3.append(ref3final)
			ref4.append(ref4final)
			ref5.append(ref5final)
			time.append(lt[0].header['JD'])

#Bin to compensate for shorter exposure times in dataset a
wasp = np.array(np.array(waspa[:len(waspa)-(len(waspa)%2)]).reshape((len(waspa)-(len(waspa)%2))/2, -1).sum(axis=1).tolist() + wasp)
ref1 = np.array(np.array(ref1a[:len(ref1a)-(len(ref1a)%2)]).reshape((len(ref1a)-(len(ref1a)%2))/2, -1).sum(axis=1).tolist() + ref1)
ref2 = np.array(np.array(ref2a[:len(ref2a)-(len(ref2a)%2)]).reshape((len(ref2a)-(len(ref2a)%2))/2, -1).sum(axis=1).tolist() + ref2)
ref3 = np.array(np.array(ref3a[:len(ref3a)-(len(ref3a)%2)]).reshape((len(ref3a)-(len(ref3a)%2))/2, -1).sum(axis=1).tolist() + ref3)
ref4 = np.array(np.array(ref4a[:len(ref4a)-(len(ref4a)%2)]).reshape((len(ref4a)-(len(ref4a)%2))/2, -1).sum(axis=1).tolist() + ref4)
ref5 = np.array(np.array(ref5a[:len(ref5a)-(len(ref5a)%2)]).reshape((len(ref5a)-(len(ref5a)%2))/2, -1).sum(axis=1).tolist() + ref5)
time = np.array(np.array(timea[:len(timea)-(len(timea)%2)]).reshape((len(timea)-(len(timea)%2))/2, -1).mean(axis=1).tolist() + time)
ref = np.mean(np.array([ref1, ref2, ref3, ref4, ref5]), axis=0)

#Calibrate to reference
cal = wasp/ref*np.median(ref)

#Note: Binning is not necessary to see the transit. This is just calculating rms deviation
nims = range(1, 16)
rms = []
nbins = []
for n in nims:
	calbin = cal[:len(cal)-(len(cal)%n)].reshape((len(cal)-(len(cal)%n))/n, -1).mean(axis=1)
	r = np.std(calbin[-int(0.1*len(calbin)):])
	nbins.append((len(cal)-(len(cal)%n))/n)
	rms.append(r)
calbinned = cal[:len(cal)-(len(cal)%9)].reshape((len(cal)-(len(cal)%9))/9, -1).mean(axis=1)
timebinned = time[:len(time)-(len(time)%9)].reshape((len(time)-(len(time)%9))/9, -1).mean(axis=1)


#Normalize to eyeballed continuum value 97000 counts because median for last 10% of points wasn't aligning at 1 for whatever reason
norm = cal/97000
normbinned = calbinned/97000

cl = []
nm = []
tm = []
for i in range(len(cal)):
	if i-5 < 0:
		b = 0
	else:
		b = i-5
	if i+5 > len(cal):
		t = len(cal)
	else:
		t = i+5
	if abs(cal[i] - np.average(cal[b:t])) < 2*np.std(cal[-int(0.1*len(cal)):]):
		cl.append(cal[i])
		nm.append(norm[i])
		tm.append(time[i])

#plt.figure()
#plt.title('WASP-43 intensity over time (uncalibrated)')
#plt.xlabel('Time (JD)')
#plt.ylabel('Intensity (counts)')
#plt.scatter(time, wasp, s=50)
#figManager = plt.get_current_fig_manager()
#figManager.window.showMaximized()
#plt.show();

#plt.figure()
#plt.title('Reference star 1 intensity over time')
#plt.xlabel('Time (JD)')
#plt.ylabel('Intensity (counts)')
#plt.scatter(time, ref1, s=50)
#figManager = plt.get_current_fig_manager()
#figManager.window.showMaximized()
#plt.show();

#plt.figure()
#plt.title('Reference star 2 intensity over time')
#plt.xlabel('Time (JD)')
#plt.ylabel('Intensity (counts)')
#plt.scatter(time, ref2, s=50)
#figManager = plt.get_current_fig_manager()
#figManager.window.showMaximized()
#plt.show();

#plt.figure()
#plt.title('Reference star 3 intensity over time')
#plt.xlabel('Time (JD)')
#plt.ylabel('Intensity (counts)')
#plt.scatter(time, ref3, s=50)
#figManager = plt.get_current_fig_manager()
#figManager.window.showMaximized()
#plt.show();

#plt.figure()
#plt.title('Reference star 4 intensity over time')
#plt.xlabel('Time (JD)')
#plt.ylabel('Intensity (counts)')
#plt.scatter(time, ref4, s=50)
#figManager = plt.get_current_fig_manager()
#figManager.window.showMaximized()
#plt.show();

#plt.figure()
#plt.title('Reference star 5 intensity over time')
#plt.xlabel('Time (JD)')
#plt.ylabel('Intensity (counts)')
#plt.scatter(time, ref5, s=50)
#figManager = plt.get_current_fig_manager()
#figManager.window.showMaximized()
#plt.show();

#plt.figure()
#plt.title('Combined reference star intensity over time')
#plt.xlabel('Time (JD)')
#plt.ylabel('Intensity (counts)')
#plt.scatter(time, ref, s=50)
#figManager = plt.get_current_fig_manager()
#figManager.window.showMaximized()
#plt.show();

#plt.figure()
#plt.title(r'WASP-43 intensity over time (bg $\sigma$-clipped, calibrated to combined reference star)')
#plt.xlabel('Time (JD)')
#plt.ylabel('Intensity (counts)')
#plt.scatter(time, cal, s=50)
#figManager = plt.get_current_fig_manager()
#figManager.window.showMaximized()
#plt.show();

#plt.figure()
#plt.title(r'WASP-43 intensity over time (bg $\sigma$-clipped, calibrated to combined reference star, 18 bins)')
#plt.xlabel('Time (JD)')
#plt.ylabel('Intensity (counts)')
#plt.scatter(timebinned, calbinned, s=50)
#figManager = plt.get_current_fig_manager()
#figManager.window.showMaximized()
#plt.show();

#plt.figure()
#plt.title(r'WASP-43 intensity over time (normalized, bg $\sigma$-clipped, calibrated to combined reference star, 18 bins)')
#plt.xlabel('Time (JD)')
#plt.ylabel('Normalized Intensity')
#plt.scatter(timebinned, normbinned, s=50)
#figManager = plt.get_current_fig_manager()
#figManager.window.showMaximized()
#plt.show();

#plt.figure()
#plt.title(r'WASP-43 intensity over time (normalized, bg $\sigma$-clipped, calibrated to combined reference star)')
#plt.xlabel('Time (JD)')
#plt.ylabel('Normalized Intensity')
#plt.scatter(time, norm, s=50)
#figManager = plt.get_current_fig_manager()
#figManager.window.showMaximized()
#plt.show();

#plt.figure()
#plt.title('Photometric precision vs number of bins for calibration to combined reference star')
#plt.xlabel('Number of bins')
#plt.ylabel('RMS deviation (counts)')
#plt.plot(nbins,rms)
#plt.gca().invert_xaxis()
#figManager = plt.get_current_fig_manager()
#figManager.window.showMaximized()
#plt.show();

#plt.figure()
#plt.title(r'WASP-43 intensity over time (bg $\sigma$-clipped, calibrated to combined reference star, data $2-\sigma$ clipped)')
#plt.xlabel('Time (JD)')
#plt.ylabel('Intensity (counts)')
#plt.scatter(tm, cl, s=50)
#figManager = plt.get_current_fig_manager()
#figManager.window.showMaximized()
#plt.show();

#plt.figure()
#plt.title(r'WASP-43 intensity over time (normalized, bg $\sigma$-clipped, calibrated to combined reference star, data $2-\sigma$ clipped)')
#plt.xlabel('Time (JD)')
#plt.ylabel('Normalized Intensity')
#plt.scatter(tm, nm, s=50)
#figManager = plt.get_current_fig_manager()
#figManager.window.showMaximized()
#plt.show();
