import matplotlib.pyplot as plt
import numpy as np

ha = 6562.8518
shift = 0.547
snr = np.array([30, 50, 70])

wv = np.arange(ha-120, ha+120, 0.31)

whwhm = 30
gsigma = 0.5
gpeaks = np.array([ha-shift, ha+shift])

wide = -15/(np.pi*whwhm*(1+((wv-ha)/whwhm)**2))
p1 = -0.1/np.sqrt(2*np.pi*gsigma**2)*np.exp(-(wv-gpeaks[0])**2/(2*gsigma**2))
p2 = -0.1/np.sqrt(2*np.pi*gsigma**2)*np.exp(-(wv-gpeaks[1])**2/(2*gsigma**2))

data = wide + p1 + p2
err = np.sqrt(abs(data)*1.45)/1.45
noise = np.random.normal(loc=data, scale=1/10000., size=len(data))
data += noise

plt.figure()
plt.title("Artificial spectrum for median separation of WD0037 and max of HE0225 (~50 km/s)")
plt.ylabel("Flux (counts)")
plt.xlabel("Wavelength (A)")
plt.plot(wv, data)
plt.xlim(6540, 6580)
plt.show()
