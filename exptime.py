import numpy as np

###Note: Goodman pxscl = 0.15, ProEM = 0.116

print("***************************\
EXPOSURE TIME CALCULATOR\
***************************")
m = float(input("Object Magnitude = "))
d = float(input("Telescope Diameter (m) = "))
pxscl = float(input("CCD Pixel Scale (\"/pix) = "))
npix = float(input("# of pixels covered by source = "))
amp = float(input("Amplitude of variability (%) = "))

mag0F = 5.e5 # # of photons from mag0 star (photons/s/cm^2)
sky = 20. #sky brightness (magnitudes/arcsecond^2)

area = np.pi * ((d*100.)/2.)**2

def bg(time, npx, pscl):
	return mag0F * 10**(-0.4 * sky) * (npx * pscl**2.) * time
def src(mag, time):
	return mag0F * 10**(-0.4 * mag) * area * time
def extime(mag, npx, pscl):
	return 2*np.sqrt(2)*(3*np.sqrt(src(mag, 1.) + bg(1., npx, pscl))/(src(mag, 1.)*amp/100.))**2.

print()
print()
print("Source photons/sec: ", src(m, 1.))
print("Variability photons/sec: ", src(m, 1.)*amp/100.)
print("Background photons/sec: ", bg(1., npix, pxscl))
print()
print("Time to reach SNR = 3: ", extime(m, npix, pxscl), "s")
