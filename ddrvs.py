###	Program to calculate RV's for DD Systems###
###	Joshua Reding UNC-CH 6/20/17	###

import numpy as np
import matplotlib.pyplot as plt
from lmfit.models import GaussianModel, LorentzianModel
from scipy.optimize import leastsq
import pyfits
import glob
#########
def systype():
	sys = int(raw_input("Single- or double-lined system? [1 or 2]: "))
	if sys < 1 or sys > 2:
		print "Unrecognized input"
		systype()
	else:
		return sys
#########
sys = systype()
spectra = glob.glob('spec*.fits')
w = 6562.8518
n = 0
rv1a = np.zeros(len(spectra))
rv1b = np.zeros(len(spectra))
rv2a = np.zeros(len(spectra))
rv2b = np.zeros(len(spectra))
for s in spectra:
	specdata = pyfits.open(s)
	fluxarr = specdata[0].data[1]
	warr = specdata[0].data[0]
	#Set the error as photon noise
	fluxerr = np.sqrt(abs(fluxarr)*1.4)/1.4
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
	cont = np.poly1d(np.polyfit(xcont, ycont, 3))
	y = np.zeros(feathigh-featlow)
	for i in range(len(y)):
		y[i] = fluxarr[featlow+i]-cont(warr[featlow+i])
	wide = LorentzianModel(prefix='wide_')
	pars = wide.guess(y, x=x)
	pars['wide_amplitude'].set(-200, max=-100)
	if sys == 1:
		rvs = np.zeros(len(spectra))
		l1 = GaussianModel(prefix='l1_')
		pars.update(l1.make_params())
		pars['l1_center'].set(w, min=w-10, max=w+10)
		pars['l1_sigma'].set(0.5, min=.1)
		pars['l1_amplitude'].set(-100, max=0)
		mod = wide+l1
		plt.figure()
		plt.plot(x, y, label='data')
		out = mod.fit(y, pars, x=x)
		comps = out.eval_components(x=x)
		print(out.fit_report(min_correl=0.5))
		plt.plot(x, out.best_fit, 'r-', label='best fit')
		plt.plot(x, comps['l1_'], 'b--', label='line1')
		plt.plot(x, comps['wide_'], 'k--', label='wide')
		plt.legend()
		plt.show()
		if w == labwavelength[0]:
			rv1a[n] = (out.params['l1_center']-w)/w*3e5
		if w == labwavelength[1]:
			rv1b[n] = (out.params['l1_center']-w)/w*3e5
	if sys == 2:
		rvs = np.zeros((2, len(spectra)))
		l1 = GaussianModel(prefix='l1_')
		pars.update(l1.make_params())
		pars['l1_center'].set(w-5, min=w-10, max=w+1)
		pars['l1_sigma'].set(0.75, min=0.5, max=1.25)
		pars['l1_amplitude'].set(-100, max=0)
		l2 = GaussianModel(prefix='l2_')
		pars.update(l2.make_params())
		pars['l2_center'].set(w+5, min=w-1, max=w+10)
		pars['l2_sigma'].set(0.25, min=0.1, max=0.75)
		pars['l2_amplitude'].set(-100, max=0)
		mod = wide+l1+l2
		plt.figure()
		plt.plot(x, y, label='data')
		out = mod.fit(y, pars, x=x)
		comps = out.eval_components(x=x)
		print(out.fit_report(min_correl=0.5))
		plt.plot(x, out.best_fit, 'r-', label='best fit')
		plt.plot(x, comps['l1_'], 'b--', label='line1')
		plt.plot(x, comps['l2_'], 'b--', label='line2')
		plt.plot(x, comps['wide_'], 'k--', label='wide')
		plt.legend()
		plt.show()
		rv1a[n] = (out.params['l1_center']-w)/w*3e5
		rv2a[n] = (out.params['l2_center']-w)/w*3e5
	n += 1

	date = specdata[0].header['DATE-OBS'].split('T')[0]
	y, m, d = date.split('-')
	time = specdata[0].header['DATE-OBS'].split('T')[1]
	h, m, s = time.split(':')


if sys == 1:
	for r in range(len(spectra)):
		rvs[r] = np.mean([rv1a[r], rv1b[r]])
if sys == 2:
	for r in range(len(spectra)):
		rvs[0, r] = np.mean([rv1a[r], rv1b[r]])
		rvs[1, r] = np.mean([rv2a[r], rv2b[r]])
print "RV1 = ", rvs[0, :]
print "RV2 = ", rvs[1, :]
