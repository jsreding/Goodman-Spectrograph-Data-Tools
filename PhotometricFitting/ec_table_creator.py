# -*- coding: utf-8 -*-
"""
Created on Fri Apr 28 13:05:18 2017
@author: Erik
"""

from ec_vizier_queries import get_wise,get_apass,get_2mass,get_galex,get_spec,get_vhs,get_atlas
import numpy as np
#import matplotlib.pyplot as plt
from astropy.table import Table 
import time

ec_ra0,ec_dec0,ec_Vmag0,BmVmag0,UmBmag0,ecname,sptype,comments=np.genfromtxt('EC_allzones_DA_decimal_allinfo.txt',delimiter=',',dtype='str',unpack=True,autostrip=True)
ec_ra,ec_dec,ec_Vmag,BmVmag,UmBmag,ecname0,sptype0,comments0=np.genfromtxt('EC_allzones_DA_decimal_allinfo.txt',delimiter=',',unpack=True,autostrip=True)
ec_Bmag = ec_Vmag + BmVmag
ec_Umag = ec_Bmag + UmBmag
    
ppmxl_dist,ppmxl_ra,ppmxl_dec,ppmxl_pmra,ppmxl_pmdec,ppmxl_b1mag,ppmxl_b2mag = np.genfromtxt('pmmxl_info.txt',dtype=float,delimiter=',',unpack=True)

all_info = []
#all_info = list(np.load("all_info.npy"))
for i in range(len(ec_ra)):
    starttime = time.time()
    print 'Now querying %s %s/%s'%(ecname[i],i+1,len(ec_ra))
    ec_eVmag = 0.05
    ec_eBmag = 0.05
    ec_eUmag = 0.2
    wise_ra0 = ppmxl_ra[i] + ((ppmxl_pmra[i]*10.5)/1000.)/(np.cos(np.deg2rad(ppmxl_dec[i]))*3600.)
    wise_dec0 = ppmxl_dec[i] + ((ppmxl_pmdec[i]*10.5)/1000.)/3600.
    apass_ra0 = ppmxl_ra[i] + ((ppmxl_pmra[i]*11.0)/1000.)/(np.cos(np.deg2rad(ppmxl_dec[i]))*3600.)
    apass_dec0 = ppmxl_dec[i] + ((ppmxl_pmdec[i]*11.0)/1000.)/3600.
    twomass_ra0 = ppmxl_ra[i] + ((ppmxl_pmra[i]*(-0.5))/1000.)/(np.cos(np.deg2rad(ppmxl_dec[i]))*3600.)
    twomass_dec0 = ppmxl_dec[i] + ((ppmxl_pmdec[i]*(-0.5))/1000.)/3600.
    galex_ra0 = ppmxl_ra[i] + ((ppmxl_pmra[i]*5.0)/1000.)/(np.cos(np.deg2rad(ppmxl_dec[i]))*3600.)
    galex_dec0 = ppmxl_dec[i] + ((ppmxl_pmdec[i]*5.0)/1000.)/3600.
    wise_sep,wise_ra,wise_dec,w1mag,ew1mag,w2mag,ew2mag = get_wise(wise_ra0,wise_dec0,2.0,local=False,ind=i)
    ap_sep,ap_ra,ap_dec,ap_B,ap_eB,ap_V,ap_eV,ap_g,ap_eg,ap_r,ap_er,ap_i,ap_ei = get_apass(apass_ra0,apass_dec0,2.0,local=False)
    tm_sep,tm_ra,tm_dec,tm_j,tm_ej,tm_h,tm_eh,tm_k,tm_ek = get_2mass(twomass_ra0,twomass_dec0,2.0,local=False,ind=i)
    gal_sep,gal_ra,gal_dec,gal_fuv,gal_efuv,gal_nuv,gal_enuv = get_galex(galex_ra0,galex_dec0,2.0,local=False)
    gia_teff,gia_logg,gia_flag,koe_teff,koe_logg,koe_flag = get_spec(ppmxl_ra[i],ppmxl_dec[i])
    at_dist,at_ra,at_dec,at_u,at_eu,at_g,at_eg,at_r,at_er,at_i,at_ei,at_z,at_ez = get_atlas(i)
    vhs_dist,vhs_ra,vhs_dec,vhs_j,vhs_ej,vhs_h,vhs_eh,vhs_k,vhs_ek = get_vhs(i)
    
# GALEX corrections: http://adsabs.harvard.edu/abs/2014MNRAS.438.3111C
    if gal_fuv and gal_nuv:    
        gal_fuv = 5.371 + (20.00*gal_fuv - 210.20)**(0.5)
        gal_nuv = 2.634 + (26.316*gal_nuv - 245.329)**(0.5)    
# ATLAS corrections: https://arxiv.org/pdf/1502.05432.pdf
    if at_dist:
        at_u = at_u - 0.01*(at_u-at_g) - 0.27
        at_g = at_g - 0.05*(at_g-at_r) + 0.06
        at_r = at_r - 0.03*(at_r-at_i) + 0.035
        at_i = at_i + 0.025
        at_z = at_z + 0.04*(at_i-at_z) - 0.04
# VISTA VHS corrections: http://casu.ast.cam.ac.uk/surveys-projects/vista/technical/photometric-properties
    if vhs_dist:
        vhs_color = vhs_j - vhs_k
        vhs_j = vhs_j + 0.031*(vhs_color)    
        vhs_h = vhs_h - 0.015*(vhs_color)
        vhs_k = vhs_k + 0.006*(vhs_color)

    all_info.append([ecname[i],sptype[i],ec_ra[i],ec_dec[i],ppmxl_dist[i],ppmxl_ra[i], \
    ppmxl_dec[i],ppmxl_pmra[i],ppmxl_pmdec[i],ppmxl_b1mag[i],ppmxl_b2mag[i],gal_sep,gal_ra,\
    gal_dec,ap_sep,ap_ra,ap_dec,at_dist,at_ra,at_dec,vhs_dist,vhs_ra,vhs_dec,tm_sep,tm_ra,tm_dec,wise_sep,wise_ra,wise_dec,gal_fuv,gal_efuv, \
    gal_nuv,gal_enuv,ec_Umag[i],ec_eUmag,ec_Bmag[i],ec_eBmag,ec_Vmag[i],ec_eVmag,ap_B,ap_eB,ap_V,ap_eV,\
    at_u,at_eu,at_g,at_eg,at_r,at_er,at_i,at_ei,at_z,at_ez,ap_g,ap_eg,ap_r,ap_er,ap_i,ap_ei, \
    vhs_j,vhs_ej,vhs_h,vhs_eh,vhs_k,vhs_ek,tm_j,tm_ej,tm_h,tm_eh,tm_k,tm_ek,w1mag,ew1mag, \
    w2mag,ew2mag,gia_teff,gia_logg,gia_flag,koe_teff,koe_logg,koe_flag])
    endtime=time.time()
    totaltime = endtime-starttime
    print 'Query took %s seconds'%totaltime

ec_table = Table(rows=all_info,names=('ec name','sptype','ec ra','ec dec','ppxml dist','ppmxl ra', \
                                'ppmxl dec','ppmxl pmra','ppmxl pmdec','ppmxl b1','ppmxl b2', \
                                'galex dist','galex ra','galex dec','apass dist','apass ra','apass dec', \
                                'atlas dist','atlas ra','atlas dec',' vhs dist','vhs ra','vhs dec',\
                                '2mass dist','2mass ra','2mass dec','wise dist','wise ra',\
                                'wise dec','galex fuv','galex efuv','galex nuv','galex enuv', \
                                'ec umag','ec eumag','ec bmag','ec ebmag','ec vmag','ec evmag', \
                                'apass bmag','apass ebmag','apass vmag','apass evmag','atlas umag',\
                                'atlas eumag','atlas gmag','atlas egmag','atlas rmag','atlas ermag',\
                                'atlas imag','atlas eimag','atlas zmag','atlas ezmag','apass gmag',\
                                'apass egmag','apass rmag','apass ermag','apass imag','apaas eimag',\
                                'vhs jmag','vhs ejmag','vhs hmag','vhs ehmag','vhs kmag','vhs ekmag',\
                                '2mass jmag','2mass ejmag','2mass hmag','2mass ehmag','2mass kmag',\
                                '2mass ekmag','wise w1mag','wise ew1mag','wise w2mag','wise ew2mag',\
                                'gia teff','gia logg','gia flag','koe teff','koe logg','koe flag'))

ec_table.write('EC_Tablev7s2r0.csv',format='csv')