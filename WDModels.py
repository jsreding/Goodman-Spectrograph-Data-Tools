import numpy as np
import csv
import sys

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]

teff = []
logg = []
M = []
mbol = []
bc = []
u = []
g = []
r = []
age = []
ug = []
gr = []
i = 0
with open(sys.argv[1], 'rb') as f:
    reader = csv.DictReader(f, delimiter=' ')
    for row in reader:
        teff.append(float(row['Teff']))
        logg.append(float(row['log g']))
        M.append(float(row['M/Mo']))
        mbol.append(float(row['Mbol']))
        bc.append(float(row['BC']))
        u.append(float(row['u']))
        g.append(float(row['g']))
        r.append(float(row['r']))
        age.append(float(row['Age']))
        ug.append(float(row['u']) - float(row['g']))
        gr.append(float(row['g']) - float(row['r']))

# print 'Teff =', teff
# print 'log g =', logg
# print 'M/Mo =', M
# print 'Mbol =', mbol
# print 'BC =', bc
# print 'u =', u
# print 'g =', g
# print 'r =', r
# print 'Age =', age
