import numpy as np
import csv
import sys

teff = []
logg = []
M = []
mbol = []
bc = []
u = []
g = []
r = []
age = []
i = 0
with open(sys.argv[1]) as f:
    reader = csv.DictReader(f)
    for row in reader:
        print row
        teff.append(row[''Teff'])
        logg.append(row['log g'])
        M.append(row['M/Mo'])
        mbol.append(row['Mbol'])
        bc.append(row['BC'])
        u.append(row['u'])
        g.append(row['g'])
        r.append(row['r'])
        age.append(row['Age'])

print 'Teff =', teff
print 'log g =', logg
print 'M/Mo =', M
print 'Mbol =', mbol
print 'BC =', bc
print 'u =', u
print 'g =', g
print 'r =', r
print 'Age =', age
