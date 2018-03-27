# -*- coding: utf-8 -*-
"""
Created on Wed Apr 26 14:58:01 2017

@author: Erik
"""

import numpy as np
from astroquery.vizier import Vizier
import astropy.coordinates as coord
import astropy.units as u

Vizier.ROW_LIMIT=-1
#Vizier.VIZIER_SERVER = 'vizier.cfa.harvard.edu'

allwise_catalog = Vizier(columns=["*","+_r"],catalog='II/328/allwise')
apass_catalog = Vizier(columns=["*","+_r"],catalog='II/336/apass9')
tmass_catalog = Vizier(columns=["*","+_r"],catalog='II/246/out')
galex_ais_catalog = Vizier(columns=["*","+_r"],catalog='II/312/ais')
galex_mis_catalog = Vizier(columns=["*","+_r"],catalog='II/312/mis')

gia_spec_catalog = Vizier(catalog='J/ApJ/743/138/table5')
koe_spec_catalog = Vizier(catalog='J/A+A/505/441/wd')

vhs_ec_cat = np.genfromtxt('vhs_query_2017-05-03.csv',delimiter=',',dtype=float)
atlas_ec_cat = np.genfromtxt('atlas_query_2017-05-03.csv',delimiter=',',dtype=float)
allwise_local = np.genfromtxt('allwise_query_2017-05-01.tbl',dtype=float)
twomass_local = np.genfromtxt('twomass_query_2017-05-01.tbl',dtype=float)
galex_local = np.genfromtxt('galex_query_2017-05-01.tbl',dtype=float,delimiter=';')
apass_local = np.genfromtxt('apass_query_2017-05-01.tbl',dtype=float,delimiter=';')

def get_wise(ra,dec,search_rad,ind=0,local=False):
    # Returns ALLWISE detection info given coords in degrees
    if local:
        cntr,wise_sep,pang,rain,decin,wise_ra,wise_dec,w1mag,ew1mag,w2mag,ew2mag = allwise_local[ind]
    else:
        search_coords = coord.SkyCoord(ra=ra,dec=dec,unit=(u.deg,u.deg),frame='icrs')
        result = allwise_catalog.query_region(search_coords,radius=coord.Angle(search_rad,"arcsec"))
        if len(result):
            wise_sep,wise_ra,wise_dec,w1mag,ew1mag,w2mag,ew2mag = result[0]['_r'][0],result[0]['RAJ2000'][0], \
                result[0]['DEJ2000'][0],result[0]['W1mag'][0],result[0]['e_W1mag'][0], \
                result[0]['W2mag'][0],result[0]['e_W2mag'][0]
        else:
            wise_sep,wise_ra,wise_dec,w1mag,ew1mag,w2mag,ew2mag = 0,0,0,0,0,0,0
    return float(wise_sep),float(wise_ra),float(wise_dec),float(w1mag),float(ew1mag),float(w2mag),float(ew2mag)

def get_apass(ra,dec,search_rad,local=False):
    # Returns ALLWISE detection info given coords in degrees
    if local:
        for i in range(len(apass_local)):
            if np.abs((ra - apass_local[i][0])) + np.abs((dec - apass_local[i][1])) < 2.77e-4:
                ra0,dec0,ap_sep,ap_ra,ap_dec,ap_V,ap_eV,ap_B,ap_eB,ap_g,ap_eg,ap_r,ap_er,ap_g,ap_eg,ap_i,ap_ei = apass_local[i]
                break
    else:
        search_coords = coord.SkyCoord(ra=ra,dec=dec,unit=(u.deg,u.deg),frame='icrs')
        result = apass_catalog.query_region(search_coords,radius=coord.Angle(search_rad,"arcsec"))
        if len(result):
            ap_sep,ap_ra,ap_dec,ap_B,ap_eB,ap_V,ap_eV,ap_g,ap_eg,ap_r,ap_er,ap_i,ap_ei = result[0]['_r'][0],result[0]['RAJ2000'][0], \
                result[0]['DEJ2000'][0],result[0]['Bmag'][0],result[0]['e_Bmag'][0], \
                result[0]['Vmag'][0],result[0]['e_Vmag'][0],result[0]['g_mag'][0], \
                result[0]['e_g_mag'][0],result[0]['r_mag'][0],result[0]['e_r_mag'][0], \
                result[0]['i_mag'][0],result[0]['e_i_mag'][0]
        else:
            ap_sep,ap_ra,ap_dec,ap_B,ap_eB,ap_V,ap_eV,ap_g,ap_eg,ap_r,ap_er,ap_i,ap_ei = 0,0,0,0,0,0,0,0,0,0,0,0,0
    return float(ap_sep),float(ap_ra),float(ap_dec),float(ap_B),float(ap_eB),float(ap_V),float(ap_eV),float(ap_g),float(ap_eg),float(ap_r),float(ap_er),float(ap_i),float(ap_ei)
    
def get_2mass(ra,dec,search_rad,ind=0,local=False):
    # Returns ALLWISE detection info given coords in degrees
    if local:
        cntr,tm_sep,pang,rain,decin,tm_ra,tm_dec,tm_j,tm_ej,jcmsig,tm_h,tm_eh,hmsig,tm_k,tm_ek,kmsig,rdflag,jh,hk,jk = twomass_local[ind]
    else:
        search_coords = coord.SkyCoord(ra=ra,dec=dec,unit=(u.deg,u.deg),frame='icrs')
        result = tmass_catalog.query_region(search_coords,radius=coord.Angle(search_rad,"arcsec"))
        if len(result):
            tm_sep,tm_ra,tm_dec,tm_j,tm_ej,tm_h,tm_eh,tm_k,tm_ek = result[0]['_r'][0],result[0]['RAJ2000'][0], \
                result[0]['DEJ2000'][0],result[0]['Jmag'][0],result[0]['e_Jmag'][0], \
                result[0]['Hmag'][0],result[0]['e_Hmag'][0],result[0]['Kmag'][0], \
                result[0]['e_Kmag'][0]
        else:
            tm_sep,tm_ra,tm_dec,tm_j,tm_ej,tm_h,tm_eh,tm_k,tm_ek = 0,0,0,0,0,0,0,0,0
    return float(tm_sep),float(tm_ra),float(tm_dec),float(tm_j),float(tm_ej),float(tm_h),float(tm_eh),float(tm_k),float(tm_ek)

def get_galex(ra,dec,search_rad,local=False):
    # Returns ALLWISE detection info given coords in degrees
    if local:
        for i in range(len(galex_local)):
            if np.abs((ra - galex_local[i][0])) + np.abs((dec - galex_local[i][1])) < 2.77e-4:
                ra0,dec0,gal_sep,gal_ra,gal_dec,gal_fuv,gal_efuv,gal_nuv,gal_enuv = galex_local[i]
                break        
    else:
        search_coords = coord.SkyCoord(ra=ra,dec=dec,unit=(u.deg,u.deg),frame='icrs')
        resulta = galex_ais_catalog.query_region(search_coords,radius=coord.Angle(search_rad,"arcsec"))
        resultb = galex_mis_catalog.query_region(search_coords,radius=coord.Angle(search_rad,"arcsec"))
        if len(resulta):
            gal_sep,gal_ra,gal_dec,gal_fuv,gal_efuv,gal_nuv,gal_enuv = resulta[0]['_r'][0],resulta[0]['RAJ2000'][0], \
                resulta[0]['DEJ2000'][0],resulta[0]['FUV'][0],resulta[0]['e_FUV'][0], \
                resulta[0]['NUV'][0],resulta[0]['e_NUV'][0]
        elif len(resultb):
            gal_sep,gal_ra,gal_dec,gal_fuv,gal_efuv,gal_nuv,gal_enuv = resultb[0]['_r'][0],resultb[0]['RAJ2000'][0], \
                resultb[0]['DEJ2000'][0],resultb[0]['FUV'][0],resultb[0]['e_FUV'][0], \
                resultb[0]['NUV'][0],resultb[0]['e_NUV'][0]
        else:
            gal_sep,gal_ra,gal_dec,gal_fuv,gal_efuv,gal_nuv,gal_enuv = 0,0,0,0,0,0,0
    return float(gal_sep),float(gal_ra),float(gal_dec),float(gal_fuv),float(gal_efuv),float(gal_nuv),float(gal_enuv)

def get_spec(ra,dec):
    search_coords = coord.SkyCoord(ra=ra,dec=dec,unit=(u.deg,u.deg),frame='icrs')
    gia_result = gia_spec_catalog.query_region(search_coords,radius=coord.Angle(10,"arcsec"))
    koe_result = koe_spec_catalog.query_region(search_coords,radius=coord.Angle(10,"arcsec"))
    gia_teff = 0
    gia_logg = 0
    koe_teff = 0
    koe_logg = 0
    gia_flag=''
    koe_flag=''
    if len(gia_result):
        gia_teff = gia_result[0]['Teff'][0]
        gia_logg = gia_result[0]['logg'][0]
        gia_flag = gia_result[0]['Notes'][0]
    if len(koe_result):
        koe_teff = koe_result[0]['Teff'][0]
        koe_logg = koe_result[0]['logg'][0]
        koe_flag = koe_result[0]['T'][0]
    return float(gia_teff),float(gia_logg),gia_flag,float(koe_teff),float(koe_logg),koe_flag
    
def get_atlas(ind):
    upload_ID,upload_RA,upload_Dec,distance,sourceid,ra,dec,uAperMag3,uAperMag3err,gAperMag3,gAperMag3err,rAperMag3,rAperMag3err,iAperMag3,iAperMag3err,zAperMag3,zAperMag3err = atlas_ec_cat[ind]
    return distance,ra,dec,uAperMag3,uAperMag3err,gAperMag3,gAperMag3err,rAperMag3,rAperMag3err,iAperMag3,iAperMag3err,zAperMag3,zAperMag3err

def get_vhs(ind):
    upload_ID,upload_RA,upload_Dec,distance,sourceID,ra,dec,jAperMag3,jAperMag3Err,hAperMag3,hAperMag3Err,ksAperMag3,ksAperMag3Err = vhs_ec_cat[ind]
    return distance,ra,dec,jAperMag3,jAperMag3Err,hAperMag3,hAperMag3Err,ksAperMag3,ksAperMag3Err