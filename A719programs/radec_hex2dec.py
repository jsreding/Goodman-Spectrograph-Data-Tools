import numpy as np
import sys

h_ra = sys.argv[1].split(':')
h_dec = sys.argv[2].split(':')

ra_h = float(h_ra[0])
ra_m = float(h_ra[1])
ra_s = float(h_ra[2])

dec_h = float(h_dec[0])
dec_m = float(h_dec[1])
dec_s = float(h_dec[2])

if ra_h >= 24 or ra_h < 0 or ra_m >= 60 or ra_m < 0 or ra_s >= 60 or ra_s < 0:
	print 'Error: Invalid RA'
	exit()

d_ra = 15 * (ra_h + ra_m / 60 + ra_s / 3600)

if dec_h > 90 or dec_h < -90 or dec_m >= 60 or dec_m < 0 or dec_s >= 60 or dec_s < 0:
	print 'Error: Invalid Dec'
	exit()

if dec_h < 0:
	dec_m = -1 * dec_m
	dec_s = -1 * dec_s

d_dec = dec_h + dec_m / 60 + dec_s / 3600

print sys.argv[1], ' = ', d_ra
print sys.argv[2], ' = ', d_dec
