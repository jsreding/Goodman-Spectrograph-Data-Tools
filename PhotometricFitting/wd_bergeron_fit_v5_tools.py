# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 11:52:26 2017

@author: Erik
"""

import numpy as np
import scipy.interpolate as interp
import scipy.optimize as opt

all_wavs = np.array([0.1528,0.2271,0.3745,0.4481,0.5423,.3540,.4470,.6222,.7632,.9049,1.235,1.662,2.159,3.352,4.602])        
all_zpts =  np.array([3631.,3631.,1649.,4060.,3723.,3631.,3631.,3631.,3631.,3631.,1594.,1024.,666.8,309.54,171.787])   
mod_wavs = np.array([all_wavs[0],all_wavs[1],all_wavs[2],all_wavs[3],all_wavs[4],all_wavs[3],all_wavs[4],all_wavs[5],\
                    all_wavs[6],all_wavs[7],all_wavs[8],all_wavs[9],all_wavs[6],all_wavs[7],all_wavs[8],all_wavs[10],\
                    all_wavs[11],all_wavs[12],all_wavs[10],all_wavs[11],all_wavs[12],all_wavs[13],all_wavs[14]])

def bergeron_DA_specsoln(logg_in,Teff_in,fit_wavs,ec_fit_flux,ec_fit_eflux,fit_bool_arr):
    Teff,logg,mass,mbol,bol_corr,jon_U,jon_B,jon_V,jon_R,jon_I,jon_J,jon_H,jon_K, \
    sdss_u,sdss_g,sdss_r,sdss_i,sdss_z,sdss_y,col_by,col_ub,col_vy,W1,W2,W3,W4 \
    ,age = np.genfromtxt('bergeron_Table_DA_WISE.txt',skip_header=4,unpack=True)    
    Teff_gal,logg_gal,mass_gal,mbol_gal,jon_V_gal,FUV,NUV,age_gal = np.genfromtxt('bergeron_Table_DA_GALEX.txt',skip_header=2,unpack=True)        
    all_data = [Teff,logg,mass,mbol,bol_corr,age,FUV,NUV,jon_U+0.0915,jon_B+0.0069,jon_V,sdss_u+0.0424,sdss_g-0.0023, \
    sdss_r-.0032,sdss_i-0.016,sdss_z-0.0276,jon_J-0.014,jon_H+0.0060,jon_K+0.008,W1,W2,W3,W4]
    Teffi = Teff[0:60]
    loggi = np.arange(7.0,9.6,0.5)    
    interp_data = []    
    for data in all_data:
        datai = np.array([data[0:60],data[60:120],data[120:180],data[180:240],data[240:300],data[300:360]])  
        biv_spline = interp.RectBivariateSpline(loggi,Teffi,datai,kx=3,ky=3)
        interp_data.append(biv_spline(logg_in,Teff_in)[0][0])     
    mmags = interp_data[6::]
    teff = Teff_in
    logg = logg_in
    mod_mags = np.array([mmags[0],mmags[1],mmags[2],mmags[3],mmags[4],mmags[3],mmags[4],mmags[5],\
                        mmags[6],mmags[7],mmags[8],mmags[9],mmags[6],mmags[7],mmags[8],mmags[10],\
                        mmags[11],mmags[12],mmags[10],mmags[11],mmags[12],mmags[13],mmags[14]])
    mod_wavs = np.array([all_wavs[0],all_wavs[1],all_wavs[2],all_wavs[3],all_wavs[4],all_wavs[3],all_wavs[4],all_wavs[5],\
                        all_wavs[6],all_wavs[7],all_wavs[8],all_wavs[9],all_wavs[6],all_wavs[7],all_wavs[8],all_wavs[10],\
                        all_wavs[11],all_wavs[12],all_wavs[10],all_wavs[11],all_wavs[12],all_wavs[13],all_wavs[14]])
    mod_zpts = np.array([all_zpts[0],all_zpts[1],all_zpts[2],all_zpts[3],all_zpts[4],all_zpts[3],all_zpts[4],all_zpts[5],\
                        all_zpts[6],all_zpts[7],all_zpts[8],all_zpts[9],all_zpts[6],all_zpts[7],all_zpts[8],all_zpts[10],\
                        all_zpts[11],all_zpts[12],all_zpts[10],all_zpts[11],all_zpts[12],all_zpts[13],all_zpts[14]])
    mod_flux = (10.**(-0.4*(mod_mags)))*mod_zpts

    mod_fit_flux = mod_flux[~fit_bool_arr]
    init_scale_ind = np.where((fit_wavs < 0.7) & (fit_wavs > 0.4))[0] 
    if len(init_scale_ind)==0:
        init_scale_ind = np.where(fit_wavs > 0)[0]
    med_scale = np.median(ec_fit_flux[init_scale_ind]/mod_fit_flux[init_scale_ind])
    init_scale_sig = np.abs(ec_fit_flux - mod_fit_flux*med_scale)/ec_fit_eflux
#    print init_scale_sig

    scale_inds = np.where((init_scale_sig < 5.0) & (fit_wavs < 3.35))[0]
#    print fit_wavs[scale_inds]
    if len(scale_inds)==0:
        scale_inds = init_scale_ind
    def fscale(wav,scl):
        return mod_fit_flux[scale_inds]*scl
    try:
        scale = opt.curve_fit(fscale,fit_wavs[scale_inds],ec_fit_flux[scale_inds],[med_scale],sigma=ec_fit_eflux[scale_inds])[0][0]
#        'print using med scale'
    except:
        scale = np.median(ec_fit_flux[scale_inds]/mod_fit_flux[scale_inds])
    mod_fit_flux_scaled = mod_fit_flux*scale
    fits=[teff,logg,-1,scale,fit_wavs,ec_fit_flux,ec_fit_eflux,mod_wavs,mod_flux,mod_fit_flux_scaled]
    return fits
    
bergeron_das = np.genfromtxt('bergeron_ecxwired_logg8.txt')
#bergeron_das = np.genfromtxt('bergeron_ecxwired.txt')

# This calculation defines a single pass bergeron grid chi-sqr calculation on the logg=8.0 grid
def fit_bergeron_grid(fit_wavs,ec_fit_flux,ec_fit_eflux,fit_bool_arr,wavlims):
    fits=[]
    chi_sqr_arr = np.zeros(len(bergeron_das))          
    ijk=0
    for mod in bergeron_das:
    # We use this block to construct flux and wav arrays of models that matches the order and length of the ec_table
        teff = mod[0]
        logg = mod[1]
        mmags = mod[6::]
        mod_mags = np.array([mmags[0],mmags[1],mmags[2],mmags[3],mmags[4],mmags[3],mmags[4],mmags[5],\
                            mmags[6],mmags[7],mmags[8],mmags[9],mmags[6],mmags[7],mmags[8],mmags[10],\
                            mmags[11],mmags[12],mmags[10],mmags[11],mmags[12],mmags[13],mmags[14]])
        mod_wavs = np.array([all_wavs[0],all_wavs[1],all_wavs[2],all_wavs[3],all_wavs[4],all_wavs[3],all_wavs[4],all_wavs[5],\
                            all_wavs[6],all_wavs[7],all_wavs[8],all_wavs[9],all_wavs[6],all_wavs[7],all_wavs[8],all_wavs[10],\
                            all_wavs[11],all_wavs[12],all_wavs[10],all_wavs[11],all_wavs[12],all_wavs[13],all_wavs[14]])
        mod_zpts = np.array([all_zpts[0],all_zpts[1],all_zpts[2],all_zpts[3],all_zpts[4],all_zpts[3],all_zpts[4],all_zpts[5],\
                            all_zpts[6],all_zpts[7],all_zpts[8],all_zpts[9],all_zpts[6],all_zpts[7],all_zpts[8],all_zpts[10],\
                            all_zpts[11],all_zpts[12],all_zpts[10],all_zpts[11],all_zpts[12],all_zpts[13],all_zpts[14]])
        mod_flux = (10.**(-0.4*(mod_mags)))*mod_zpts

        teff = mod[0]
        logg = mod[1]
        mod_fit_flux = mod_flux[~fit_bool_arr]
        chi_inds = np.where((fit_wavs < wavlims[1]) & (fit_wavs > wavlims[0]))[0]
        init_scale_ind = np.where((fit_wavs < 0.7) & (fit_wavs > 0.4))[0] 
        if len(init_scale_ind)==0:
            init_scale_ind = np.where(fit_wavs > 0)[0]
        med_scale = np.median(ec_fit_flux[init_scale_ind]/mod_fit_flux[init_scale_ind])
        init_scale_sig = np.abs(ec_fit_flux - mod_fit_flux*med_scale)/ec_fit_eflux
#        print init_scale_sig
        scale_inds = np.where((init_scale_sig < 5.0) & (fit_wavs < 3.35))[0]
#        print fit_wavs[scale_inds]
        if len(scale_inds)==0:
            scale_inds = init_scale_ind
        def fscale(wav,scl):
            return mod_fit_flux[scale_inds]*scl
        try:
            scale = opt.curve_fit(fscale,fit_wavs[scale_inds],ec_fit_flux[scale_inds],[med_scale],sigma=ec_fit_eflux[scale_inds])[0][0]
            'print using med scale'
        except:
            scale = np.median(ec_fit_flux[scale_inds]/mod_fit_flux[scale_inds])
        mod_fit_flux_scaled = mod_fit_flux*scale
        chi_sqr = np.sum(((mod_fit_flux_scaled[chi_inds] - ec_fit_flux[chi_inds])**2.0)/(ec_fit_eflux[chi_inds]**2.))
        chi_sqr_arr[ijk] = chi_sqr/(len(ec_fit_flux[chi_inds])-1)
        fits.append([teff,logg,chi_sqr/(len(ec_fit_flux[chi_inds])-1),scale,fit_wavs,ec_fit_flux,ec_fit_eflux,mod_wavs,mod_flux,mod_fit_flux_scaled])
        ijk+=1
    if len(np.where(chi_sqr_arr==chi_sqr_arr.min())[0]):
        best_fit = np.where(chi_sqr_arr==chi_sqr_arr.min())[0][0]
    else:
        best_fit=0
    return fits[best_fit]
