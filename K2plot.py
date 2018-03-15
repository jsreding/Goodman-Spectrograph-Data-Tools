import sys
import matplotlib.pyplot as plt
import numpy as np
import nufft
import glob

f = open(sys.argv[1], "r")
lines = f.readlines()
objname = lines[0].split(' ')[1]+lines[0].split(' ')[2]
time = np.zeros(len(lines)-1)
flux = np.zeros(len(lines)-1)
if 'go' in sys.argv[1] == True:
    unc = np.zeros(len(lines)-1)
i = 0
for l in lines[1:]:
    data = l.split(' ')
    time[i] = float(data[0])
    flux[i] = float(data[1])*100
    if 'go' in sys.argv[1] == True:
        unc[i] = float(data[2])
    i += 1

freq = np.linspace(2*np.pi*.1e-6, 2*np.pi*310.e-6, 10000)
ft = abs(nufft.nufft1d3(time, flux, freq))
fq = freq[ft[35:].argmax()+35]/(2*np.pi)
per = 1/fq/86400
# print "Highest peak:", bestper/86400, "days, ", freq[mx]*1000000, "microHertz"
# phase = np.fmod(time, bestper)

fig = plt.figure()
ax1 = fig.add_subplot(311)
ax1.minorticks_on()
ax1.set_title(objname)
ax1.set_ylabel('Rel. Flux (%)')
ax1.set_xlabel('Time (days)')
ax1.set_xlim([0, 80])
if 'go' in sys.argv[1] == True:
    ax1.errorbar(time/86400., flux, yerr=unc, fmt='o', linestyle="None", ms=2)
else:
    ax1.scatter(time/86400., flux)
ax2 = fig.add_subplot(312)
ax2.minorticks_on()
ax2.text(0.65, 0.9,r'Highest Peak: $%s days (%s ppt), %s \mu Hz$'%(np.round(per, 3), np.round(ft[35:].max()*20, 3), np.round(fq*1e6, 3)), horizontalalignment='center', verticalalignment='center', transform=ax2.transAxes, fontsize=14)
ax2.set_ylabel('Amplitude (ppt)')
ax2.set_xlabel(r'Possible frequencies ($\mu Hz$)')
ax2.set_xlim([0, 310])
ax2.plot(freq*1000000/(2*np.pi), ft*20)
ax3 = fig.add_subplot(313)
ax3.minorticks_on()
ax3.set_ylabel('Amplitude (ppt)')
ax3.set_xlabel(r'Possible frequencies ($\mu Hz$)')
ax3.set_xlim([0, 50])
ax3.plot(freq*1000000/(2*np.pi), ft*20)
plt.show()
