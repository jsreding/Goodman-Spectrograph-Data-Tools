import numpy as np
import csv
import sys
from scipy import interpolate as interp
import matplotlib.pyplot as plt
import time

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
        M.append(float(row['M/Mo']))
        mbol.append(float(row['Mbol']))
        bc.append(float(row['BC']))
        u.append(float(row['u']))
        g.append(float(row['g']))
        r.append(float(row['r']))
        age.append(float(row['Age']))
        ug.append((float(row['u'])+0.0424) - (float(row['g'])-0.0023))
        gr.append((float(row['g'])-0.0023) - (float(row['r'])-0.0032))

def colorfit(teff_in, logg_in):
    teffi = teff[:58]
    loggi = np.arange(7.0, 9.6, 0.5)
    color = [ug, gr]
    intdata = []
    for data in color:
        datai = np.array([data[:58], data[58:116], data[116:174], data[174:232], data[232:290], data[290:348]])
        bvspl_teff = interp.RectBivariateSpline(loggi, teffi, datai)
        intdata.append(bvspl_teff(logg_in, teff_in)[0][0])
    return intdata

tr = np.linspace(1500, 120000, 1000)
lgr = np.linspace(7.0, 9.5, 1000)
TR, LGR = np.meshgrid(tr, lgr)
UG = colorfit(TR, LGR)[0]
print UG.shape()

def simplex(ugi, gri):
    n = 2
    ugi = np.float(ugi)
    gri = np.float(gri)
    col = np.array([ugi, gri])
    for c in range(len(col)):
        t = []
        l = []
        t_mid = 65000
        t_diff = 50000
        l_mid = 8.25
        l_diff = 1.25
        e = 1
        if c == 0:
            print "Calculating curve for u-g"
        elif c == 1:
            print "Calculating curve for g-r"
        while len(t) < 100:
            t_test = np.random.uniform(t_mid-t_diff, t_mid+t_diff)
            l_test = np.random.uniform(l_mid-l_diff, l_mid+l_diff)
            test = colorfit(t_test, l_test)[c]
            if abs(test - col[c]) < 10**(-e):
                t_mid = t_test
                t_diff = t_diff/10.
                l_mid = l_test
                l_diff = l_diff/2.
                e += 1
                if e == 4:
                    t.append(t_test)
                    l.append(l_test)
                    t_mid = 65000
                    t_diff = 50000
                    l_mid = 8.25
                    l_diff = 1.25
                    e = 1
        if c == 0:
            ugt = t
            ugl = l
        elif c == 1:
            grt = t
            grl = l

    ugfunc = np.polyfit(ugt, ugl, 3)
    grfunc = np.polyfit(grt, grl, 3)

    Tfinal = np.roots(ugfunc-grfunc)[0]
    return Tfinal, np.poly1d(ugfunc)(Tfinal)

    # t = np.linspace(1500, 120000, 1000)
    #
    # plt.figure()
    # plt.xlabel("Teff")
    # plt.ylabel("log g")
    # plt.xlim(1500, 120000)
    # plt.ylim(7., 9.5)
    # plt.plot(t, np.poly1d(ugfunc)(t), label='u-g =%s'%(ugi))
    # plt.plot(t, np.poly1d(grfunc)(t), label='g-r =%s'%(gri))
    # plt.legend()
    # plt.show()

# input = 'n'
# while input != 'y':
#     T = raw_input("Teff: ")
#     L = raw_input("log g: ")
#     print "u-g:", colorfit(T, L)[0]
#     print "g-r:", colorfit(T, L)[1]
#     input = raw_input("Done? (y/n) ")


# trange = np.linspace(1500, 120000, 1000)
# lgrange = np.arange(7.0, 9.6, 0.5)
# plt.figure()
# plt.title("Teff, log g surface in u-g, g-r parameter space")
# plt.ylabel("u-g")
# plt.xlabel("g-r")
# plt.xlim(-0.6, 0.5)
# plt.ylim(-0.6, 1)
# for l in lgrange:
#     ugr = []
#     grr = []
#     for t in trange:
#         ugr.append(colorfit(t, l)[0])
#         grr.append(colorfit(t, l)[1])
#     plt.plot(grr, ugr, label="l=%s"%l, zorder=1)
# tcoarse = np.arange(1500, 120100, 500)
# for t in tcoarse:
#     ugr = []
#     grr = []
#     for l in lgrange:
#         ugr.append(colorfit(t, l)[0])
#         grr.append(colorfit(t, l)[1])
#     plt.plot(grr, ugr, color='Black', zorder=1)
# plt.gca().invert_yaxis()
# plt.show()

# dahteff = [14078., 15254., 17326., 20335., 11140., 16750., 11351., 22775., 19134., 17904., 17936., 22490., 25636., 8753., 6886., 16175., 18764., 18529., 15613., 15980., 10604., 18090., 22510., 12284., 10182., 19142., 18694., 9871., 7525., 20160., 13458., 12275., 14592., 16700., 28315., 47547., 66229., 28125., 21696., 59983., 61801., 10620.]
# dahlogg = [8.45, 9.02, 9.02, 9.0, 8.738, 9.12, 8.33, 9.05, 7.81, 7.79, 8.21, 8.69, 7.44, 9.06, 9.0, 8.37, 9.16, 8.09, 8.15, 8.16, 7.97, 7.91, 10.0, 8.4, 6.84, 7.86, 8.07, 8.04, 9.02, 9.0, 9.0, 9.0, 5.65, 8.32, 10.0, 8.84, 7.58, 9.22, 8.08, 8.0, 6.799, 9.189]
# ugr = []
# grr = []
# for i in range(len(dahteff)):
#     ugr.append(colorfit(dahteff[i], dahlogg[i])[0])
#     grr.append(colorfit(dahteff[i], dahlogg[i])[1])
# plt.scatter(grr, ugr, label='DAH', marker='o', zorder=2)
#
# dbahteff = [18254., 19992., 17772., 18421., 16005., 23629., 15413.]
# dbahlogg = [8.07, 8.02, 8.17, 7.85, 8.13, 8.32, 8.02]
# ugr = []
# grr = []
# for i in range(len(dbahteff)):
#     ugr.append(colorfit(dbahteff[i], dbahlogg[i])[0])
#     grr.append(colorfit(dbahteff[i], dbahlogg[i])[1])
# plt.scatter(grr, ugr, label='DBAH', marker='^', zorder=2)
#
# dbhteff = [34521., 40000.]
# dbhlogg = [8.44, 9.61]
# ugr = []
# grr = []
# for i in range(len(dbhteff)):
#     ugr.append(colorfit(dbhteff[i], dbhlogg[i])[0])
#     grr.append(colorfit(dbhteff[i], dbhlogg[i])[1])
# plt.scatter(grr, ugr, label='DBH', marker='s', zorder=2)
#
# dbzhmteff = [16503.]
# dbzhmlogg = [8.33]
# ugr = []
# grr = []
# for i in range(len(dbzhmteff)):
#     ugr.append(colorfit(dbzhmteff[i], dbzhmlogg[i])[0])
#     grr.append(colorfit(dbzhmteff[i], dbzhmlogg[i])[1])
# plt.scatter(grr, ugr, label='DBZHM', marker='x', zorder=2)
#
# dqhteff = [11700.]
# dqhlogg = [8.5]
# ugr = []
# grr = []
# for i in range(len(dqhteff)):
#     ugr.append(colorfit(dqhteff[i], dqhlogg[i])[0])
#     grr.append(colorfit(dqhteff[i], dqhlogg[i])[1])
# plt.scatter(grr, ugr, label='DQH', marker='*', zorder=2)
# plt.xlim(-0.6, 0.5)
# plt.ylim(1.0, -0.6)
# plt.legend()
# plt.show()

# t1 = raw_input("u-g ")
# t2 = raw_input("g-r ")
# T = np.zeros(100)
# L = np.zeros(100)
# tm = np.zeros(100)
# for s in range(100):
#     start = time.time()
#     T[s], L[s] = simplex(t1, t2)
#     print "--- %s seconds ---"%(time.time() - start)
#     tm[s] = time.time() - start
#
# print np.mean(tm)
# print np.sum(tm)
# print np.mean(T), "+/-", np.std(T)
# print np.mean(L), "+/-", np.std(L)
# plt.figure()
# plt.title("Scatter in 100 trials of MC technique")
# plt.xlabel("Teff")
# plt.ylabel("log g")
# plt.scatter(T, L)
# plt.show()
