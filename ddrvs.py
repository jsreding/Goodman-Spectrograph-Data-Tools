###	Program to calculate RV's for DD Systems###
###		Joshua Reding UNC-CH 6/20/17		###

import numpy as np
import matplotlib.pyplot as plt
from lmfit.models import GaussianModel, LorentzianModel, PolynomialModel
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
	contlow = np.argmax(warr > (w - 150.))
	conthigh = np.argmax(warr > (w - 100.))
	featlow = np.argmax(warr > (w - 20.))
	feathigh = np.argmax(warr > (w + 20.))
	#Fit Lorentizan to wide feature and Gaussian(s) to small
	x = warr[featlow:feathigh]
	ycont = fluxarr[contlow:conthigh]
	xcont = warr[contlow:conthigh]
	contmod = PolynomialModel(3, prefix='cont_')
	contpars = contmod.guess(ycont, x=xcont)
	cont = np.mean(contmod.fit(ycont, contpars, x=xcont).eval(x=x))
	y = fluxarr[featlow:feathigh]-cont
	wide = LorentzianModel(prefix='wide_')
	pars = wide.guess(y, x=x)
	if sys == 1:
		rvs = np.zeros(len(spectra))
		l1 = GaussianModel(prefix='l1_')
		pars.update(l1.make_params())
		pars['l1_center'].set(w, min=w-10, max=w+10)
		pars['l1_sigma'].set(0.5, min=.1)
		pars['l1_amplitude'].set(-300, max=0)
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
		pars['l1_center'].set(np.random.normal(w, 1)-2.5, min=w-10, max=w-0.5)
		pars['l1_sigma'].set(0.5, min=0.1, max=1)
		pars['l1_amplitude'].set(-300, max=0)
		l2 = GaussianModel(prefix='l2_')
		pars.update(l2.make_params())
		pars['l2_center'].set(np.random.normal(w, 1)+2.5, min=w+0.5, max=w+10)
		pars['l2_sigma'].set(0.5, min=0.1, max=1)
		pars['l2_amplitude'].set(-300, max=0)
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
		# if w == labwavelength[0]:
		rv1a[n] = (out.params['l1_center']-w)/w*3e5
		rv2a[n] = (out.params['l2_center']-w)/w*3e5
		# if w == labwavelength[1]:
		# 	rv1b[n] = (out.params['l1_center']-w)/w*3e5
		# 	rv2b[n] = (out.params['l2_center']-w)/w*3e5
	n += 1

if sys == 1:
	for r in range(len(spectra)):
		rvs[r] = np.mean([rv1a[r], rv1b[r]])
if sys == 2:
	for r in range(len(spectra)):
		rvs[0, r] = np.mean([rv1a[r], rv1b[r]])
		rvs[1, r] = np.mean([rv2a[r], rv2b[r]])
print rvs
