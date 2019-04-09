import matplotlib.pyplot as plt
import numpy as np
from astropy.convolution import convolve, Box1DKernel

from matplotlib import rc, rcParams
rcParams['font.family'] = 'afm'
rcParams['font.sans-serif'] = ['Helvetica']
rc('text', usetex=True)

f = open("./SDSSJ1252-0234toSDSSg.dat")
r = open("./SDSSJ1252_redtoSDSSr.dat")
mod = open("./Koester8500_8")
lines = f.readlines()
rlines = r.readlines()
mlines = mod.readlines()
wave = np.zeros(len(lines))
flux = np.zeros(len(lines))
rwave = np.zeros(len(rlines))
rflux = np.zeros(len(rlines))
mwave = np.zeros(len(mlines))
mflux = np.zeros(len(mlines))
n = 0
for l in lines:
    if not l.startswith('#'):
        wave[n] = float(l.strip().split(' ')[0])
        flux[n] = float(l.strip().split(' ')[1])
        n += 1
n = 0
for l in rlines:
    if not l.startswith('#'):
        rwave[n] = float(l.strip().split(' ')[0])
        rflux[n] = float(l.strip().split(' ')[1])
        n += 1
n = 0
for l in mlines:
    if not l.startswith('#'):
        mwave[n] = float(l.strip().split(' ')[0])
        mflux[n] = float(l.strip().split(' ')[1])/8.75e22
        n += 1

wave = np.hstack((wave, rwave))
flux = np.hstack((flux, rflux))

flux = convolve(flux, Box1DKernel(7))

fig = plt.figure()
ax = fig.add_subplot(111)
ax.spines['top'].set_color('none')
ax.spines['bottom'].set_color('none')
ax.spines['left'].set_color('none')
ax.spines['right'].set_color('none')
ax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
ax.set_ylabel(r"$f_\lambda$ (erg/cm$^2$/s/\AA)", fontsize=20)
ax.set_xlabel(r"Wavelength (\AA)", fontsize=20)
ax1 = fig.add_subplot(211)
ax1.set_ylim(1.5e-16, 6.5e-16)
ax1.set_xlim(3675., 5300)
ax1.plot(mwave, mflux, color='C1')
ax1.plot(wave, flux, color='C0')
ax2 = fig.add_subplot(212)
ax2.set_xlim(5675., 7240.)
ax2.set_ylim(1.e-16, 4.e-16)
ax2.plot(mwave, mflux, color='C1', label=r"Koester Model; $\mathrm{T_{eff}}=8500, \log g = 8.0$")
ax2.plot(wave, flux, color='C0', label="Observed")
ax2.legend(fontsize=14)
plt.show()
