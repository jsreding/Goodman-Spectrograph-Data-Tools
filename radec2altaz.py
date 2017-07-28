#Created by Joshua Reding for ASTR 719, Spring 2017
#Converts RA and Dec of object to Alt and Az for given location, date, and time
#Call in format 'python radec2altaz.py year month day time(HH:MM:SS) obs_long obs_lat RA(HH:MM:SS) Dec(DD:MM:SS)'

#Call relevant packages
import numpy as np
import datetime as dt
import sys

#Read input arguments
y = int(sys.argv[1])
m = int(sys.argv[2])
d = int(sys.argv[3])
time = sys.argv[4].split(':')
lng = float(sys.argv[5])
lat = float(sys.argv[6])
h_ra = sys.argv[7].split(':')
h_dec = sys.argv[8].split(':')

#Separate hex values
ra_h = float(h_ra[0])
ra_m = float(h_ra[1])
ra_s = float(h_ra[2])

dec_h = float(h_dec[0])
dec_m = float(h_dec[1])
dec_s = float(h_dec[2])

time_h = int(time[0])
time_m = int(time[1])
time_s = int(time[2])

#Convert hex RA and Dec to decimal degrees, time to decimal hour
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

if time_h >= 24 or time_h < 0 or time_m >= 60 or time_m < 0 or time_s >= 60 or time_s < 0:
	print 'Error: Invalid Time'
	exit()
d_time = time_h + time_m / 60.0 + time_s / 3600.0

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
t_now = dt.datetime(year=y, month=m, day=d, hour=time_h, minute=time_m, second=time_s)
days = t_now - J2000
d_days = days.days + float(days.seconds/(24.0*60.0*60.0))
lst = 100.46 + 0.985647*d_days + lng + 15.0*d_time

#Bring LST into range 0-360
while lst < 0:
	lst += 360
while lst > 360:
	lst -= 360

#Calculate HA, bring into range 0-360
ha = lst - d_ra
while ha < 0:
	ha += 360
while ha > 360:
	ha -= 360

#Calculate altitude
alt = np.arcsin(np.sin(np.radians(d_dec))*np.sin(np.radians(lat))+np.cos(np.radians(d_dec))*np.cos(np.radians(lat))*np.cos(np.radians(ha)))*180/np.pi
print 'alt=', np.round(alt, 3)

#Calculate azimuth
a = np.arccos((np.sin(np.radians(d_dec))-np.sin(np.radians(alt))*np.sin(np.radians(lat)))/(np.cos(np.radians(alt))*np.cos(np.radians(lat))))*180/np.pi
if np.sin(np.radians(ha)) >= 0:
	az = 360 - a
if np.sin(np.radians(ha)) < 0:
	az = a
print 'az=', np.round(az, 3)
