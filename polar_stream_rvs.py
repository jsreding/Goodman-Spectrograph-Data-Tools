'''
Program written September 2015 by JTF.

This program attempts to fit multiple gaussians to emission lines. Written primarilly for fitting polar accretion streams.

'''

import numpy as np
import matplotlib.pyplot as plt
import mpfit
import pyfits as fits
import os
from scipy.special import erf

#Define fitting functions
def gauss_skew(x,p):
    return p[3] + (2./p[0]) * p[4] * np.exp(-((x-p[1])/p[0])**2./2.) * (1. + erf((p[2]*(x-p[1])/p[0])/np.sqrt(2.)))/2. + p[4]*np.exp(-(((x-p[5])/(np.sqrt(2)*p[6])))**2.)
#    p[0] = sigma
#    p[1] = x_naught
#    p[2] = skew
#    x = data

def fitskew(p,fjac=None,x=None,y=None,err=None):
    model = gauss_skew(x,p)
    status = 0
    return([status,(y-model)/err])

def gauss(x,p): #This is if you only want to plot a single gaussian
    return p[0] + p[1]*x +  p[2]*np.exp(-(((x-p[3])/(np.sqrt(2)*p[4])))**2.)
# 

def fitgauss(p,fjac=None,x=None,y=None,err=None):
    #Parameter values are passed in p
    #fjac = None just means partial derivatives will not be computed
    model = gauss(x,p)
    status = 0
    return([status,(y-model)/err])

def multigauss(x,p):
    return p[0] + p[1]*x + p[2]*np.exp(-(((x-p[3])/(np.sqrt(2)*p[4])))**2.) + p[5]*np.exp(-(((x-p[6])/(np.sqrt(2)*p[7])))**2.)

def fitmultigauss(p,fjac=None,x=None,y=None,err=None):
    #Parameter values are passed in p
    #fjac = None just means partial derivatives will not be computed
    model = multigauss(x,p)
    status = 0
    return([status,(y-model)/err])

#Change directory to where spectra are
os.chdir('/afs/cas.unc.edu/depts/physics_astronomy/soar/pipeline/2014/2014-05-03/GOODMAN/J1928/400line/')

#Lab wavelength of the line you want to fit
labwavelength = 6678.151 #6296.5, 6562.79, 6678.151, 4685.804

#Read in file with list of spectral names.
spectra = np.genfromtxt('list_stream.txt',dtype=str)

result = np.zeros([len(spectra),5])
resulterr = np.zeros([len(spectra),5])
skyresult = np.zeros([len(spectra),8])
skyresulterr = np.zeros([len(spectra),8])
offset = np.zeros(len(spectra))
velocity = np.zeros(len(spectra))
velocityerr = np.zeros(len(spectra))

percent = np.zeros(len(spectra))

n = 0.
for spec in spectra:
    #Read in spectrum
    print spec
    dataval = fits.open(spec)
    fluxarr = dataval[0].data[0,0,:]

    #Set the error as photon noise
    fluxerr = np.sqrt(fluxarr*1.4)/1.4

    #Set up wavelengths
    specwav0 = dataval[0].header['CRVAL1']
    specdeltawav = dataval[0].header['CD1_1']
    warr = np.zeros(len(fluxarr))
    warr[0] = specwav0
    ival = np.arange(1,len(fluxarr))
    for i in ival:
        warr[i] = warr[i-1] + specdeltawav

    #if on first iteration, use spectra to create guesses for mpfit. Otherwise, use previous result
    if n >= 0:
        #guess is the array of guess for the fitting
        indexlow = np.argmax(warr > (labwavelength - 30.))
        indexmed = np.argmax(warr > (labwavelength - 20.))
        indexhigh = np.argmax(warr > (labwavelength + 10.))
        guess = np.zeros(5)
        guess[0] = np.mean(fluxarr[indexlow:indexmed]) #this is the continuum value
        guess[1] = (fluxarr[indexmed]-fluxarr[indexlow])/(warr[indexmed]-warr[indexlow])#slope of continuum
        guess[2] = np.max(fluxarr[indexmed:indexhigh]) - guess[0] #max value for broad component
        guess[3] = labwavelength #Central value for broad component
        guess[4] = 11. #sigma of broad component
        #guess[5] = guess[1] #peak of narrow component
        #guess[6] = labwavelength + 6. #central value of narrow component
        #guess[7] = 2. #sigma of central component
       # guess = np.array([5.,6563.,0.,50.,50.,6570.,2.])
    else:
        #Use the previous solved parameters as a starting point
        guess = fitparams.params

    #constrain parameters for fit
    #limited is boolean on whether or not lower/upper limits exist on a particular parameter.
    #limits sets the limits if limited = 1
    #The initial guess must be within this limit otherwise it will not work.
    params = [{'limits':[0,0],'limited':[0,0]} for i in range(5)] #8 total parameters
    #params[4]['limited'] = [1.,0.] 
    #params[3]['limited'] = [1,0]
    #params[6]['limited'] = [1,1]
    #params[4]['limits'] = [1.,0.] 
    #params[6]['limits'] = [0.5,5.] 
    #params[1]['limited'] = [0.,1.]
    #params[1]['limits'] = [0.,-1.]

    fithigh = np.argmax(warr > (labwavelength + 100.))
    fitlow = np.argmax(warr > (labwavelength - 60.))
    lambdas = warr[fitlow:fithigh]
    fluxes = fluxarr[fitlow:fithigh]
    err = fluxerr[fitlow:fithigh]
    fa = {'x':lambdas,'y':fluxes,'err':err}
    #fitparams = mpfit.mpfit(fitmultigauss,guess,functkw=fa,parinfo=params)
    fitparams = mpfit.mpfit(fitgauss,guess,functkw=fa)

    #Save the fitted parameters
    result[n,:] = fitparams.params
    resulterr[n,:] = fitparams.perror
    
    #Use this part if you want to measure/correct for a skyline
    skyarr = dataval[0].data[1,0,:]
    if n ==0:
        skyguess = [16.,10.,90.,5890.7,0.9,55.,5896.75,0.9] #Continuum, max, central value, sigma,max,central value, sigma
    #    skyguess = [40.,10.,800.,5577.,7.]
    else:
        skyguess = skyparams.params

    if spec == 'bt0473.J1928_1200_series2.ms.fits':
        skyguess = [16.,10.,90.,5890.7,0.9,55.,5896.75,0.9] #Continuum, max, central value, sigma,max,central value, sigma

    skylow = np.argmax(warr > 5850.)
    skyhigh = np.argmax(warr > 5925.)
    skylambda = warr[skylow:skyhigh]
    skyflux = skyarr[skylow:skyhigh]
    skyerr = np.sqrt(skyflux)
    sa = {'x':skylambda,'y':skyflux,'err':skyerr}
    skyparams = mpfit.mpfit(fitmultigauss,skyguess,functkw=sa)

    #print spec
    #plt.clf()
    #plt.plot(skylambda,skyflux)
    #plt.plot(skylambda,multigauss(skylambda,skyparams.params))
    #plt.show()
    
    skyresult[n,:] = skyparams.params
    skyresulterr[n,:] = skyparams.perror

    offset[n] = ((5889.953 - skyresult[n,3]) + (5895.923 - skyresult[n,6]))/2.
    #offset[n] = 5577.338 - skyresult[n,3]
    result[n,3] = result[n,3] + offset[n]
    #result[n,5] = result[n,5] + offset[n]
    

    #Convert from wavelength to velocity
    velocity[n] = (result[n,3] - labwavelength)/labwavelength * 2.99792e5
    velocityerr[n] = (resulterr[n,3])/labwavelength * 2.99792e5
    

    #Let's plot!
    #broadparams = [fitparams.params[0],fitparams.params[1],fitparams.params[2],fitparams.params[3],fitparams.params[4]]
    #broadfit = gauss(lambdas,broadparams)
    #narrowparams = [fitparams.params[0],fitparams.params[1],fitparams.params[5],fitparams.params[6],fitparams.params[7]]
    #narrowfit = gauss(lambdas,narrowparams)

    #Compute how much absorption take up
    #continuum = fitparams.params[0] + fitparams.params[1]*fitparams.params[3]
    #absorp = gauss(fitparams.params[3],fitparams.params)
    #percent[n] = (1. - absorp/continuum)*100.
    #print continuum, absorp
    #print percent
    
    print spec, offset[n]
    print fitparams.params
    print fitparams.perror
    plt.clf()
    lambdafit = np.linspace(6620.0,6778.0,num=500)
    plt.plot(lambdas,fluxes,'b')
    fit = gauss(lambdafit,fitparams.params)
    plt.plot(lambdafit,fit,'g',linewidth=2.0)
    #plt.plot(lambdas,broadfit,'r')
    #plt.plot(lambdas,narrowfit,'c')
    #plt.show()



    n += 1.

print len(velocity)
x = np.linspace(1,59,59)
plt.clf()
#plt.plot(percent)
#plt.errorbar(x,result[:,2],yerr=resulterr[:,2])
plt.errorbar(x,velocity,yerr=1.5*velocityerr)
#plt.plot(x,offset)
plt.show()

othervel, othererr, newphase = np.genfromtxt('J1928_05-03_400_Halpha_stream_vel.txt',unpack=True)

#Save results
#np.savetxt('J1928_05-03_1200_He_stream_vel.txt',np.transpose([velocity,1.5*velocityerr,newphase]))
#np.savetxt('J1928_05-03_400_abs_dip.txt',np.transpose([percent]))
