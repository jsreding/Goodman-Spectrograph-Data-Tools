###	Program to calculate RV's for DD Systems###
###	Joshua Reding UNC-CH 6/20/17	###

import numpy as np
import matplotlib.pyplot as plt
from lmfit.models import GaussianModel, LorentzianModel
import scipy.signal
from scipy.optimize import leastsq
import pyfits
import glob
from astropy.time import Time

w = 6562.8518

#########
def systype():
	sys = int(raw_input("Single- or double-lined system? [1 or 2]: "))
	if sys < 1 or sys > 2:
		print "Unrecognized input"
		systype()
	else:
		return sys
def maxrv():
	mxrv = int(raw_input("Maximum expected RV: "))
	mxshift = mxrv/3e5*w
	print "Maximum expected wavelength shift =", mxshift
	return mxshift, mxrv

def onclick(event):
	global c
	global h
	c.append(event.xdata)
	h.append(event.ydata)
	if len(c) == 2:
		plt.close()
#########

sys = systype()
mxshift, mxrv = maxrv()
spectra = glob.glob('spec*.fits')
n = 0
rv1 = []
rv2 = []
rv = []
jd = []
imnum = np.zeros(len(spectra))
for s in spectra:
	imnum[n] = int(s.split('_')[1])
	specdata = pyfits.open(s)
	fluxarr = specdata[0].data[1]
	warr = specdata[0].data[0]
	#Define continuum and feature ranges
	contlowa = np.argmax(warr > (w - 180.))
	contlowb = np.argmax(warr > (w - 130.))
	featlow = np.argmax(warr > (w - 120.))
	feathigh = np.argmax(warr > (w + 120.))
	conthigha = np.argmax(warr > (w + 130.))
	conthighb = np.argmax(warr > (w + 180.))
	#Fit Lorentizan to wide feature and Gaussian(s) to small
	x = warr[featlow:feathigh]
	ycont = np.array(np.ndarray.tolist(fluxarr[contlowa:contlowb])+np.ndarray.tolist(fluxarr[conthigha:conthighb]))
	xcont = np.array(np.ndarray.tolist(warr[contlowa:contlowb])+np.ndarray.tolist(warr[conthigha:conthighb]))
	cont = np.poly1d(np.polyfit(xcont, ycont, 2))
	y = np.zeros(feathigh-featlow)
	for i in range(len(y)):
		y[i] = (fluxarr[featlow+i]-cont(warr[featlow+i]))/cont(warr[featlow+i])
	#Set the error as photon noise, generate fake spectrum
	gain = specdata[0].header['GAIN']
	fluxerr = np.sqrt(abs(y)*gain)/gain
	fspec = np.random.normal(loc=y, scale=fluxerr, size=len(y))
	wide = LorentzianModel(prefix='wide_')
	pars = wide.guess(y, x=x)
	pars['wide_center'].set(w, min=w-10, max=w+10)
	pars['wide_amplitude'].set(-17000, max=0)
	if sys == 1:
		rvs = np.zeros(len(spectra))
		l1 = GaussianModel(prefix='l1_')
		pars.update(l1.make_params())
		pars['l1_center'].set(w, min=w-10, max=w+10)
		pars['l1_sigma'].set(0.5, min=.1)
		pars['l1_amplitude'].set(-.3, max=0)
		mod = wide+l1
		plt.figure()
		plt.plot(x, y, label='data')
		# plt.plot(x, fspec, label='fake data')
		out = mod.fit(y, pars, x=x)
		comps = out.eval_components(x=x)
		print(out.fit_report(min_correl=0.5))
		plt.plot(x, out.best_fit, 'r-', label='best fit')
		plt.plot(x, comps['l1_'], 'b--', label='line1')
		plt.plot(x, comps['wide_'], 'k--', label='wide')
		plt.legend()
		plt.show()
		rv1[n] = (out.params['l1_center']-w)/w*3e5
	if sys == 2:
		l1 = GaussianModel(prefix='l1_')
		pars.update(l1.make_params())
		c = []
		h = []
		fig = plt.figure()
		# plt.xlabel('Wavelength (A)')
		# plt.ylabel('Flux (counts)')
		plt.plot(x, y, label='data')
		# plt.plot(x, fspec, label='fakedata')
		plt.legend()
		plt.xlim(6540, 6580)
		plt.title(s)
		cid = fig.canvas.mpl_connect('button_press_event', onclick)
		plt.show()
		pars['l1_center'].set(c[0], min=c[0]-0.5, max=c[0]+0.5)
		pars['l1_sigma'].set(0.3, min=0.2, max=0.5)
		pars['l1_amplitude'].set(h[0], max=0)
		# pars['l1_height'].set(max=-200)
		l2 = GaussianModel(prefix='l2_')
		pars.update(l2.make_params())
		pars['l2_center'].set(c[1], min=c[1]-0.5, max=c[1]+0.5)
		pars['l2_sigma'].set(0.3, min=0.2, max=0.5)
		pars['l2_amplitude'].set(h[1], max=0)
		# pars['l2_height'].set(max=-200)
		mod = wide+l1+l2
		out = mod.fit(y, pars, x=x)
		comps = out.eval_components(x=x)
		print(out.fit_report(min_correl=0.5))
		plt.figure()
		plt.xlabel('Wavelength (A)')
		plt.ylabel('Flux (counts)')
		plt.plot(x, y, label='data')
		# plt.plot(x, fspec, label='fake data')
		plt.xlim(6540, 6580)
		plt.title(s)
		plt.plot(x, out.best_fit, 'r-', label='best fit')
		# plt.plot(x, comps['l1_'], 'b--', label='line1')
		# plt.plot(x, comps['l2_'], 'b--', label='line2')
		# plt.plot(x, comps['wide_'], 'k--', label='wide')
		plt.legend()
		plt.show()
		r1 = (out.params['l1_center']-w)/w*3e5
		r2 = (out.params['l2_center']-w)/w*3e5
		if abs(out.params['l1_amplitude']) > abs(out.params['l2_amplitude']):
			rv1.append(r1)
			rv2.append(r2)
			rv.append(r1-r2)
		else:
			rv1.append(r2)
			rv2.append(r1)
			rv.append(r2-r1)
		time = specdata[0].header['DATE-OBS']
		t = Time(time, format='isot', scale='utc')
		jd.append(t.jd-2450000.)
	n += 1

if sys == 1:
	for r in range(len(rv)):
		rvs[r] = rv1[r]
	print "RV1 = ", rvs[:]
	print "JD-2450000 = ", jd
if sys == 2:
	jd = np.array(jd)
	rv = np.array(rv)
	rv1 = np.array(rv1)
	rv2 = np.array(rv2)
	print "RV1 = ", rv1
	print "RV2 = ", rv2
	print "RV = ", rv
	print "JD = ", jd

	t = np.linspace(jd[0], jd[-1], 10000)
	freq = np.linspace(2*np.pi/0.1, 2*np.pi/10, 10000)
	normval = len(jd)
	pgram = scipy.signal.lombscargle(jd, rv, freq)
	mx = np.argmax(pgram)
	bestper = 2*np.pi/freq[mx]
	print bestper
	phase = np.fmod(t, bestper)
	eqn1 = np.zeros(len(t))
	eqn2 = np.zeros(len(t))
	for e in range(len(t)):
		eqn1[e] = mxrv*np.sin(freq[mx]*t[e])
		eqn2[e] = -mxrv*np.sin(freq[mx]*t[e])

	plt.figure()
	plt.title('Lomb-Scargle Periodogram')
	plt.ylabel('Power')
	plt.xlabel('Period (days)')
	plt.plot(2*np.pi/freq, np.sqrt(4*(pgram/normval)))
	plt.show()

	plt.figure()
	plt.xlabel('JD (-2450000)')
	plt.ylabel('RV (km/s)')
	plt.plot(t, eqn1, label='fit1')
	plt.plot(t, eqn2, label='fit2')
	plt.scatter(jd, rv1, label='rv1')
	plt.scatter(jd, rv2, label='rv2')
	plt.scatter(jd, rv, label=r'$\Delta RV$')
	plt.legend()
	plt.show()
