#Created by Joshua Reding for ASTR 719, Spring 2017
#Converts RA and Dec of object to Alt for given location, date, and plots altitude over time
#Call in format 'python radec2altaz.py year month day obs_lat obs_long listfile'

#Call relevant packages
import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
import sys

file = open(sys.argv[6], "r")

for line in file.readlines():

	data = line.split('	')
	objname = str(data[0])

	#Read input arguments
	y = int(sys.argv[1])
	m = int(sys.argv[2])
	d = int(sys.argv[3])
	lat = float(sys.argv[4])
	lng = float(sys.argv[5])
	h_ra = data[1].split(':')
	h_dec = data[2].split(':')

	#Separate hex values
	ra_h = float(h_ra[0])
	ra_m = float(h_ra[1])
	ra_s = float(h_ra[2])

	dec_h = float(h_dec[0])
	dec_m = float(h_dec[1])
	dec_s = float(h_dec[2])

	#Convert hex RA and Dec to decimal degrees
	if ra_h >= 24 or ra_h < 0 or ra_m >= 60 or ra_m < 0 or ra_s >= 60 or ra_s < 0:
		print 'Error: Invalid RA'
		exit()
	d_ra = 15 * (ra_h + ra_m / 60.0 + ra_s / 3600.0)

	if dec_h > 90 or dec_h < -90 or dec_m >= 60 or dec_m < 0 or dec_s >= 60 or dec_s < 0:
		print 'Error: Invalid Dec'
		exit()
	if dec_h < 0:
		dec_m = -1 * dec_m
		dec_s = -1 * dec_s
	d_dec = dec_h + dec_m / 60.0 + dec_s / 3600.0

	#Check validity of latitude and longitude, convert negative longitudes to range 0-360 if necessary
	if lat > 90.0 or lat < -90.0:
		print 'Error: Invalid Latitude'
		exit()

	if lng > 360.0 or lng < -360.0:
		print 'Error: Invalid Longitude'
		exit()
	if lng < 0:
		lng += 360

	#Calculate LST (correction values taken from http://www.stargazing.net/kepler/altaz.html)
	J2000 = dt.datetime(year=2000, month=1, day=1, hour=12, minute=0, second=0)
	t = np.linspace(0, 23.999, 100)
	time_h = np.zeros(100)
	time_m = np.zeros(100)
	time_s = np.zeros(100)
	d_days = np.zeros(100)
	lst = np.zeros(100)
	ha = np.zeros(100)
	alt = np.zeros(100)
	for i in range(0, 100):
		time_h[i] = int(t[i])
		time_m[i] = int((t[i] - time_h[i])*60)
		time_s[i] = int(((t[i] - time_h[i])*60 - time_m[i])*60)
		d_days[i] = (dt.datetime(year=y, month=m, day=d, hour=int(time_h[i]), minute=int(time_m[i]), second=int(time_s[i])) - J2000).days + float((dt.datetime(year=y, month=m, day=d, hour=int(time_h[i]), minute=int(time_m[i]), second=int(time_s[i])) - J2000).seconds/(24.0*60.0*60.0))
		lst[i] = 100.46 + 0.985647*d_days[i] + lng + 15.0*t[i]
	#Bring LST into range 0-360
		while lst[i] < 0:
			lst[i] += 360
		while lst[i] > 360:
			lst[i] -= 360
	#Calculate HA, bring into range 0-360
		ha[i] = lst[i] - d_ra
		while ha[i] < 0:
			ha[i] += 360
		while ha[i] > 360:
			ha[i] -= 360
	#Calculate altitude
		alt[i] = np.degrees(np.arcsin(np.sin(np.radians(d_dec))*np.sin(np.radians(lat))+np.cos(np.radians(d_dec))*np.cos(np.radians(lat))*np.cos(np.radians(ha[i]))))

	#Plot altitude over time
	plt.plot(t, alt, label = '%s' %(objname));

plt.xlabel('Time (UT)')
plt.xlim([0, 24])
plt.ylabel('Altitude (degrees)')
plt.ylim([0, 90])
plt.title('Altitude of objects at location Lat=%s Long=%s on date %s/%s/%s' %(sys.argv[5], sys.argv[4], sys.argv[1], sys.argv[2], sys.argv[3]))
plt.xticks(np.arange(0, 24))
plt.grid(True)
plt.plot([0,24],[30,30], 'k-', lw = 2, label = 'Airmass = 2.0');
plt.legend();
plt.show()
