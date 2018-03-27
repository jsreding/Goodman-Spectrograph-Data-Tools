# -*- coding: utf-8 -*-
"""
Created on Tue Aug 09 15:45:04 2016

@author: Erik
"""

import numpy as np
from astropy.table import Table 
import matplotlib.pyplot as plt
import os,glob
from wd_bergeron_fit_v5_tools import bergeron_DA_specsoln,fit_bergeron_grid

plt.ioff()

dir_name='test_plots'
try:
    os.makedirs(dir_name)
except:
    pass

#files = glob.glob(dir_name+'/*')
#for f in files:
#    os.remove(f)

all_wavs = np.array([0.1528,0.2271,0.3745,0.4481,0.5423,.3540,.4470,.6222,.7632,.9049,1.235,1.662,2.159,3.352,4.602])        
all_zpts =  np.array([3631.,3631.,1649.,4060.,3723.,3631.,3631.,3631.,3631.,3631.,1594.,1024.,666.8,309.54,171.787])   
mod_wavs = np.array([all_wavs[0],all_wavs[1],all_wavs[2],all_wavs[3],all_wavs[4],all_wavs[3],all_wavs[4],all_wavs[5],\
                    all_wavs[6],all_wavs[7],all_wavs[8],all_wavs[9],all_wavs[6],all_wavs[7],all_wavs[8],all_wavs[10],\
                    all_wavs[11],all_wavs[12],all_wavs[10],all_wavs[11],all_wavs[12],all_wavs[13],all_wavs[14]])

ec_table = Table.read('EC_Tablev3r0.csv')
wd_Teff_fits = np.zeros(len(ec_table))
wd_logg_fits = np.zeros(len(ec_table))
wd_Teff_phot = np.zeros(len(ec_table))
wd_logg_phot = np.zeros(len(ec_table))
wd_chisqr_fits = np.zeros(len(ec_table))
wd_scale_fits = np.zeros(len(ec_table))
excess_flag = np.zeros(len(ec_table))
ec_fluxtbl = []
ec_fluxtbl_err = []
ec_wavstbl = []
wd_mod_fluxtbl = []
wd_mod_wavstbl = []
wise_expected = np.zeros(len(ec_table))
wise_good_det = np.zeros(len(ec_table))
test_inds = [443]
for i in test_inds:#range(len(ec_table)):
    print 'Now fitting: EC '+ec_table[i]['ec name']
    ec_mags = np.array([ec_table[i]['galex fuv'],ec_table[i]['galex nuv'],ec_table[i]['ec umag'],ec_table[i]['ec bmag'],ec_table[i]['ec vmag'],ec_table[i]['apass bmag'], \
                ec_table[i]['apass vmag'],-1,-1,-1,-1,-1,ec_table[i]['apass gmag'],ec_table[i]['apass rmag'], \
                ec_table[i]['apass imag'],ec_table[i]['vhs jmag'],ec_table[i]['vhs hmag'],ec_table[i]['vhs kmag'], \
                ec_table[i]['2mass jmag'],ec_table[i]['2mass hmag'],ec_table[i]['2mass kmag'],ec_table[i]['wise w1mag'],\
                ec_table[i]['wise w2mag']])
#    ec_mags = ec_mags- np.array([0,0,0.083,0.074,0.055,0.074,0.055,0.0,0.0,0.0,0.0,0.0,0.067,0.046,0.035,0.015,0.009,0.006,0.015,0.009,0.006,0.0002,0.003])
#    ec_mags = np.array([ec_table[i]['galex fuv'],ec_table[i]['galex nuv'],ec_table[i]['ec umag'],ec_table[i]['ec bmag'],ec_table[i]['ec vmag'],ec_table[i]['apass bmag'], \
#                ec_table[i]['apass vmag'],ec_table[i]['atlas umag'],ec_table[i]['atlas gmag'],ec_table[i]['atlas rmag'],ec_table[i]['atlas imag'],ec_table[i]['atlas zmag'],ec_table[i]['apass gmag'],ec_table[i]['apass rmag'], \
#                ec_table[i]['apass imag'],ec_table[i]['vhs jmag'],ec_table[i]['vhs hmag'],ec_table[i]['vhs kmag'], \
#                ec_table[i]['2mass jmag'],ec_table[i]['2mass hmag'],ec_table[i]['2mass kmag'],ec_table[i]['wise w1mag'],\
#                ec_table[i]['wise w2mag']])
    ec_emags = np.array([ec_table[i]['galex efuv'],ec_table[i]['galex enuv'],ec_table[i]['ec eumag'],ec_table[i]['ec ebmag'],ec_table[i]['ec evmag'],ec_table[i]['apass ebmag'], \
                ec_table[i]['apass evmag'],ec_table[i]['atlas eumag'],ec_table[i]['atlas egmag'],ec_table[i]['atlas ermag'],ec_table[i]['atlas eimag'],ec_table[i]['atlas ezmag'],ec_table[i]['apass egmag'],ec_table[i]['apass ermag'], \
                ec_table[i]['apaas eimag'],ec_table[i]['vhs ejmag'],ec_table[i]['vhs ehmag'],ec_table[i]['vhs ekmag'], \
                ec_table[i]['2mass ejmag'],ec_table[i]['2mass ehmag'],ec_table[i]['2mass ekmag'],ec_table[i]['wise ew1mag'],\
                ec_table[i]['wise ew2mag']])

#  Quick filter to handle any left over weird magnitudes
    for l in range(len(ec_mags)):
        if ec_mags[l] > 20.0 or ec_mags[l] < 1.0 or ec_emags[l] > 20.0 or ec_emags[l] < 0.0:
            ec_mags[l] = np.nan
            ec_emags[l] = np.nan

# Flux and flux error propogations
    ec_zpts = np.array([all_zpts[0],all_zpts[1],all_zpts[2],all_zpts[3],all_zpts[4],all_zpts[3],all_zpts[4],all_zpts[5],\
                        all_zpts[6],all_zpts[7],all_zpts[8],all_zpts[9],all_zpts[6],all_zpts[7],all_zpts[8],all_zpts[10],\
                        all_zpts[11],all_zpts[12],all_zpts[10],all_zpts[11],all_zpts[12],all_zpts[13],all_zpts[14]])
    ec_flux = (10.**(-0.4*(ec_mags)))*ec_zpts
    ec_eflux = np.log(10)*0.4*ec_flux*ec_emags

# Now we proceed by filtering out nans in mags and emags
    fit_bool_arr = np.logical_or(np.isnan(ec_mags),np.isnan(ec_emags)) #boolean array that finds nans in the mags OR emags array in question

    fit_wavs = mod_wavs[~fit_bool_arr]
    ec_fit_mags = ec_mags[~fit_bool_arr]
    ec_fit_emags = ec_emags[~fit_bool_arr]
    ec_fit_flux = ec_flux[~fit_bool_arr]
    ec_fit_eflux = ec_eflux[~fit_bool_arr]
    
# quick check for reasonable flux errors ( > 5%)
    for l in range(len(ec_fit_flux)):
        if ec_fit_eflux[l]/ec_fit_flux[l] < 0.05:
            ec_fit_eflux[l] = ec_fit_flux[l]*0.05

    wavlims = [0.0,1.23]
# first we calculate a chi-sqr grid using all wavs blue of J with the scale anchored to optical phot. 
    best_fit = fit_bergeron_grid(fit_wavs,ec_fit_flux,ec_fit_eflux,fit_bool_arr,wavlims) 
    chi_sqr1 = best_fit[2]
#        print best_fit[0],best_fit[1],best_fit[2]
# Check to see if IR excess at WISE wavelengths
    excess_inds = np.where(best_fit[4] > 2.16)[0]
    excess_check = (best_fit[5][excess_inds] - best_fit[9][excess_inds])/best_fit[6][excess_inds]
    phot_fit = best_fit    
    if len(excess_check):    
        if any(excess_check > 5.0) or all(excess_check > 3.0):
            excess_flag[i]=1
# check spec flag to determine if we need to find a photometric solution
    spec_fit = np.zeros(len(phot_fit))
    if ec_table[i]['gia teff'] > 0:
        print 'Using spectroscopic solution from Gianninas'
        spec_fit = bergeron_DA_specsoln(ec_table[i]['gia logg'],ec_table[i]['gia teff'],fit_wavs,ec_fit_flux,ec_fit_eflux,fit_bool_arr)
        excess_inds = np.where(spec_fit[4] > 2.16)[0]
        excess_check = (spec_fit[5][excess_inds] - spec_fit[9][excess_inds])/spec_fit[6][excess_inds]
        if len(excess_check):
            if any(excess_check > 5.0) or all(excess_check > 3.0):
                excess_flag[i] = 1
    if ec_table[i]['koe teff'] > 0:
        print 'Using spectroscopic solution from Koester'
        spec_fit = bergeron_DA_specsoln(ec_table[i]['koe logg'],ec_table[i]['koe teff'],fit_wavs,ec_fit_flux,ec_fit_eflux,fit_bool_arr)
        excess_inds = np.where(spec_fit[4] > 2.16)[0]
        excess_check = (spec_fit[5][excess_inds] - spec_fit[9][excess_inds])/spec_fit[6][excess_inds]
        if len(excess_check):
            if any(excess_check > 5.0) or all(excess_check > 3.0):
                excess_flag[i] = 1    
    spec_fit[0]=0
    if spec_fit[0]:
        final_fit = spec_fit
    else:
        final_fit = phot_fit

    wd_Teff_phot[i]=phot_fit[0]
    wd_logg_phot[i]=phot_fit[1]

    wd_Teff_fits[i]=final_fit[0]
    wd_logg_fits[i]=final_fit[1]
    wd_chisqr_fits[i]=final_fit[2]
    wd_scale_fits[i]=final_fit[3]
    
    pwav = final_fit[4]
    pflux = final_fit[5]
    pflux_err = final_fit[6]
    mod_fit_flux = final_fit[9]
    all_mod_wavs = final_fit[7]
    all_mod_flux = final_fit[8]*wd_scale_fits[i]
#    if all_mod_flux[-2] > 0.000080 or all_mod_flux[-1] > 0.00011:
#        wise_expected[i] = 1
#    snr_flux = pflux/pflux_err
#    if snr_flux[-1] > 5 or snr_flux[-2] >5:
#        wise_good_det[i]=1
    sort_inds = np.argsort(all_mod_wavs)
    name = ec_table[i]['ec name']
    otype = ec_table[i]['sptype']
    if len(pwav):
        plt.figure()
        plt.clf()
        if excess_flag[i]:
            plt.title(name+' Type:%s'%otype + ' Chi Sqr:%s'%np.round(wd_chisqr_fits[i],decimals=1) + ' IR Excess!')
        else:
            plt.title(name+' Type:%s'%otype + ' Chi Sqr:%s'%np.round(wd_chisqr_fits[i],decimals=1))
        plt.xscale('log')
        plt.yscale('log')
        plt.errorbar(pwav,pflux,yerr=pflux_err,fmt='o',color='b',ecolor='b',markeredgecolor='b',markersize=2)
        plt.plot(all_mod_wavs[sort_inds],all_mod_flux[sort_inds],'r--',label='Logg %s Teff %s'%(np.round(wd_Teff_fits[i],decimals=1),np.round(wd_logg_fits[i],decimals=1)))
#        plt.plot([3.35,4.6],[0.000080,0.00011],'rv')        
        plt.xlabel('Wavelength (micron)')
        plt.xlim(0.1,11.)
        plt.ylabel('Flux (Jy)')
        plt.legend(loc=8)
        fig1 = plt.gcf()
        fname = dir_name +'/'+ name +'.jpg'
        fig1.savefig(fname)
        plt.close()

#ec_table['excess flag'] = Table.Column(excess_flag,name='excess flag')
#ec_table['wd teff'] = Table.Column(wd_Teff_fits,name='wd teff')
#ec_table['wd logg'] = Table.Column(wd_logg_fits,name='wd logg')
#ec_table['wd phot teff'] = Table.Column(wd_Teff_phot,name='wd phot teff')
#ec_table['wd phot logg'] = Table.Column(wd_logg_phot,name='wd phot logg')
#ec_table['wd chisqr'] = Table.Column(wd_chisqr_fits,name='wd chisqr')
#ec_table['wd flux scale'] = Table.Column(wd_scale_fits,name='wd flux scale')

#ec_table.write('EC_Tablev3r1.csv',format='csv')   