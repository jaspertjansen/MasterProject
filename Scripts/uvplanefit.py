# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 13:34:57 2024

@author: jansen
"""

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import RegscorePy

FWHM = '2FWHM'
Amp_on = False

#Title = 'AS2COS0008.1'; vis='{}_NRAO_combined_target.ms.split.{}'.format(Title, FWHM)
#Title = 'AS2COS0013.1'; vis='{}_NRAO_target.ms.split.{}'.format(Title, FWHM)
#Title = 'AS2COS0014.1'; vis='{}_NRAO_all_combined_target_nobckgrndsrc.ms.split.{}'.format(Title, FWHM)
#Title = 'AS2COS0028.1'; vis='{}_NRAO_target.ms.split.{}'.format(Title, FWHM); vel=640
#Title = 'AS2COS0031.1'; vis='{}_NRAO_target_nobckgrndsrc.ms.split.{}'.format(Title, FWHM)
#Title = 'AS2COS0054.1'; vis='{}_NRAO_target.ms.split.{}'.format(Title, FWHM)
#Title = 'AS2COS0066.1'; vis='{}_NRAO_combined_target.ms.split.{}'.format(Title, FWHM)
#Title = 'AS2COS0139.1'; vis='{}_NRAO_target.ms.split.{}'.format(Title, FWHM)
#Title = 'AS2UDS010.0'; vis='{}_Marta_combined_target.ms.split.{}'.format(Title, FWHM)
Title = 'AS2UDS012.0'; vis='{}_NRAO_target.ms.split.{}'.format(Title, FWHM); vel=523
#Title = 'AS2UDS126.0'; vis='{}_NRAO_combined_target.ms.split.{}'.format(Title, FWHM)
#Title = 'AEG2'; vis='{}_Marta_target.ms.split.{}'.format(Title, FWHM)
#Title = 'CDFN1'; vis='{}_NRAO_target.ms.split.{}'.format(Title, FWHM)
#Title = 'CDFN2'; vis='{}_NRAO_target.ms.split.{}'.format(Title, FWHM)
#Title = 'CDFN8'; vis='{}_NRAO_target_nobckgrndsrc.ms.split.{}'.format(Title, FWHM)


#Title = 'AS2UDS012.0'; vis='{}_NRAO_target.ms.split.{}.mock'.format(Title, FWHM)


if FWHM=='2FWHM':
    vel=vel*2


number_max_fit = 60 #1000

uvbinsize = 5

UVDataFile = 'D:\\Master Astronomy Research year 2\\Master Project\\uvfit\\'+vis+'_binned_'+(str)(uvbinsize)+'_klambda.txt'

def gauss1d_new (x, a, err):
    return a* np.exp(-(x**2)/(2*(err)**2))

def straight2(x, B): # this is your 'straight line' y=f(x)
    return 0*x+B

def convert_uv_to_lm(theta_uv): #fromhttps://casadocs.readthedocs.io/en/stable/notebooks/UVTaper_Imaging_Weights.html
    theta_lm  = 3600 * ( 4*np.log(2)/np.pi ) / (theta_uv * np.pi/180.0)
    return theta_lm

data_file = open(UVDataFile)
UVDist, Real, RealError, Imag, ImagError, NoVis = np.loadtxt(data_file, unpack = True)
data_file.close()

UVDist = UVDist[~np.isnan(UVDist)]
Real = Real [~np.isnan(Real)]
RealError = RealError [~np.isnan(RealError)]
Imag = Imag [~np.isnan(Imag)]
ImagError = ImagError [~np.isnan(ImagError)]

Amp = (Real**2+Imag**2)**(1/2)
AmpError = RealError * (Real/Amp) + ImagError * (Imag/Amp)

pred_params_real, pcov_real = curve_fit(gauss1d_new, UVDist[UVDist<number_max_fit], Real[UVDist<number_max_fit], sigma = RealError[UVDist<number_max_fit], maxfev=100000, bounds=([0, -np.inf], [0.5, np.inf]))
pred_err_real = np.sqrt(np.diag(pcov_real))
residual_real = Real - gauss1d_new(UVDist, *pred_params_real)
chisq_real = sum((residual_real[UVDist<number_max_fit] / RealError[UVDist<number_max_fit]) ** 2)
print("Chisq Real Gaus", round(chisq_real))

if Amp_on:
    pred_params_amp, pcov_amp = curve_fit(gauss1d_new, UVDist[UVDist<number_max_fit], Amp[UVDist<number_max_fit], sigma = AmpError[UVDist<number_max_fit], maxfev=100000, bounds=([0, -np.inf], [0.5, np.inf]))
    pred_err_amp = np.sqrt(np.diag(pcov_amp))
    residual_amp = Amp - gauss1d_new(UVDist, *pred_params_amp)
    chisq_amp = sum((residual_amp[UVDist<number_max_fit] / AmpError[UVDist<number_max_fit]) ** 2)
    print("Chisq Amp Gaus", round(chisq_amp))

pred_params_real2, pcov_real2 = curve_fit(straight2, UVDist[UVDist<number_max_fit], Real[UVDist<number_max_fit], sigma = RealError[UVDist<number_max_fit], maxfev=100000, bounds=([0], [0.5]))
pred_err_real2 = np.sqrt(np.diag(pcov_real2))
residual_real2 = Real - straight2(UVDist, *pred_params_real2)
chisq_real2 = sum((residual_real2[UVDist<number_max_fit] / RealError[UVDist<number_max_fit]) ** 2)
print("Chisq Real point source", round(chisq_real2))

if Amp_on:
    pred_params_amp2, pcov_amp2 = curve_fit(straight2, UVDist[UVDist<number_max_fit], Amp[UVDist<number_max_fit], sigma = AmpError[UVDist<number_max_fit], maxfev=100000, bounds=([0], [0.5]))
    pred_err_amp2 = np.sqrt(np.diag(pcov_amp2))
    residual_amp2 = Amp - straight2(UVDist, *pred_params_amp2)
    chisq_amp2 = sum((residual_amp2[UVDist<number_max_fit] / AmpError[UVDist<number_max_fit]) ** 2)
    print("Chisq Amp point source", round(chisq_amp2))

#pred_err_real = [10000,10000]

print ("Real R_1/2 =",round(pred_params_real[1]), '+/-', round(pred_err_real[1]), 'klambda')

arcsec_real = convert_uv_to_lm(pred_params_real[1]*1e3)
arcsec_err_real = (pred_err_real[1]/pred_params_real[1])*convert_uv_to_lm(pred_params_real[1]*1e3)
print ("Real R_1/2 =", round(arcsec_real,1), '+/-', round(arcsec_err_real,1), 'arcsec')

print ("Real Gaus total flux =",round(pred_params_real[0],2), '+/-', round(pred_err_real[0],2), 'mJy')

print ("Real point source total flux=",round(pred_params_real2[0],2), '+/-', round(pred_err_real2[0],2), 'mJy')

if Amp_on:
    print ("Amp R_1/2 =",round(pred_params_amp[1]), '+/-', round(pred_err_amp[1]), 'klambda')
    
    arcsec_amp = convert_uv_to_lm(pred_params_amp[1]*1e3)
    arcsec_err_amp = (pred_err_amp[1]/pred_params_amp[1])*convert_uv_to_lm(pred_params_amp[1]*1e3)
    print ("Amp R_1/2 =", round(arcsec_amp,1), '+/-', round(arcsec_err_amp,1), 'arcsec')
    
    print ("Amp Gaus total flux =",round(pred_params_amp[0],2), '+/-', round(pred_err_amp[0],2), 'mJy')
    
    print ("Amp point source total flux=",round(pred_params_amp2[0],2), '+/-', round(pred_err_amp2[0],2), 'mJy')


aic_gaus = RegscorePy.aic.aic(Real[UVDist<number_max_fit], gauss1d_new(UVDist, pred_params_real[0], pred_params_real[1])[UVDist<number_max_fit],2)
aic_straight = RegscorePy.aic.aic(Real[UVDist<number_max_fit], straight2(UVDist, pred_params_real2)[UVDist<number_max_fit],1)

bic_gaus = RegscorePy.bic.bic(Real[UVDist<number_max_fit], gauss1d_new(UVDist, pred_params_real[0], pred_params_real[1])[UVDist<number_max_fit],2)
bic_straight = RegscorePy.bic.bic(Real[UVDist<number_max_fit], straight2(UVDist, pred_params_real2)[UVDist<number_max_fit],1)

plot_UVDist = np.linspace(0, max(UVDist),10000)

plt.errorbar(UVDist, Real, yerr=RealError, fmt='.', color='black', elinewidth=1)

plt.plot(plot_UVDist, gauss1d_new(plot_UVDist, pred_params_real[0], pred_params_real[1]), color='blue', linestyle='dashed', \
         label='chisq={}, {} +- {} mJy, {} +- {} "'.format(round(chisq_real), round(pred_params_real[0],2), round(pred_err_real[0],2), round(arcsec_real,1),round(arcsec_err_real,1)))

plt.plot(plot_UVDist, straight2(plot_UVDist, pred_params_real2), color='blue', label='chisq={}, {} +- {} mJy'.format(round(chisq_real2), round(pred_params_real2[0],2), round(pred_err_real2[0],2)))

if Amp_on:
    plt.errorbar(UVDist, Amp, yerr=AmpError, fmt='.', color='orange')
    plt.plot(plot_UVDist, gauss1d_new(plot_UVDist, pred_params_amp[0], pred_params_amp[1]), color='orange', linestyle='dashed', \
             label='chisq={}, {} +- {} mJy, {} +- {} "'.format(round(chisq_amp), round(pred_params_amp[0],2), round(pred_err_amp[0],2), round(arcsec_amp,1),round(arcsec_err_amp,1)))
    plt.plot(UVDist, straight2(UVDist, pred_params_amp2), color='orange', label='chisq={}, {} +- {} mJy'.format(round(chisq_amp2), round(pred_params_amp2[0],2), round(pred_err_amp2[0],2)))

plt.plot([0,110],[0,0], c = '0', lw = 1, linestyle = 'dashed', zorder = 3)

plt.xlabel('UVdist (klambda)')
plt.ylabel('Real (mJy)')
#plt.ylim(bottom=-0.1)
plt.text(5,0.7, Title+", $\Delta v={}$ km/s".format(vel), verticalalignment='center', fontsize = 15)

if aic_gaus<aic_straight and bic_gaus<bic_straight:

    plt.text(5,0.60, "Gaussian profile", verticalalignment='center', fontsize = 10)
    plt.text(5,0.55, "$I_{{max}}=({}\\pm{})$ mJy".format(round(pred_params_real[0],2), round(pred_err_real[0],2)), verticalalignment='center', fontsize = 10)
    plt.text(5,0.50, "$R_{{1/2}}=({}\\pm{})$ arscec".format(round(arcsec_real,1), round(arcsec_err_real,1)), verticalalignment='center', fontsize = 10)

else:
    plt.text(5,0.60, "Point source", verticalalignment='center', fontsize = 10)
    plt.text(5,0.55, "$I_{{max}}=({}\\pm{})$ mJy".format(round(pred_params_real2[0],2), round(pred_err_real2[0],2)), verticalalignment='center', fontsize = 10)


plt.ylim(-0.25,0.75)
plt.xlim(0,110)
#plt.legend()
#plt.title('{}, Ignoring data after 60klambda for fit and chisq\n Blue=Real, Orange=Amp'.format(Title))
#plt.savefig('D:\\Master Astronomy Research year 2\\Master Project\\uvplanefit_{}_{}.png'.format(FWHM, Title))
plt.show()















##%%
#
#
#import numpy as np
#from scipy.optimize import curve_fit
#
#import matplotlib.pyplot as plt
#
#vis = 'AS2UDS012.0_NRAO_target.ms.split.line.shifted'
##vis = 'AS2COS0028.1_NRAO_target.ms.split.line.shifted'
##vis = 'AS2COS0013.1_NRAO_target.ms.split.line'
##vis = 'AS2COS0014.1_NRAO_all_combined_target.ms.split.line'
#
#UVDataFile = 'D:\\Master Astronomy Research year 2\\Master Project\\uvfit\\'+vis+'_binned_varbinwidth_klambda.txt'
#
#def gauss1d_new (x, a, err):
#    return a* np.exp(-(x**2)/(2*(err)**2))
#
#
#def convert_uv_to_lm(theta_uv): #fromhttps://casadocs.readthedocs.io/en/stable/notebooks/UVTaper_Imaging_Weights.html
#    theta_lm  = 3600 * ( 4*np.log(2)/np.pi ) / (theta_uv * np.pi/180.0)
#    return theta_lm
#
#data_file = open(UVDataFile)
#UVDist, Real, RealError, Imag, ImagError, NoVis = np.loadtxt(data_file, unpack = True)
#data_file.close()
#
##NoVis = NoVis[~np.isnan(UVDist)]
#UVDist = UVDist[~np.isnan(UVDist)]
#Real = Real[~np.isnan(Real)]
#RealError = RealError[~np.isnan(RealError)]
#Imag = Imag[~np.isnan(Imag)]
#ImagError = ImagError[~np.isnan(ImagError)]
#
#Amp = (Real**2+Imag**2)**(1/2)
#AmpError = RealError * (Real/Amp) + ImagError * (Imag/Amp)
#
#
#pred_params_real, pcov_real = curve_fit(gauss1d_new, UVDist, Real, sigma = RealError, p0= [0.15,30], maxfev=100000)
#pred_params_real, pcov_real = curve_fit(gauss1d_new, UVDist, Real, sigma = RealError, p0= [0.15,4], maxfev=100000)
#pred_params_amp, pcov_amp = curve_fit(gauss1d_new, UVDist, Amp, sigma = AmpError, p0= [0.3,10], maxfev=100000)
#
#
#pred_err_real = np.sqrt(np.diag(pcov_real))
#pred_err_amp = np.sqrt(np.diag(pcov_amp))
#
#
##print UVDataFile,
#print ("Real R_1/2 =",pred_params_real[1], '+/-', pred_err_real[1], 'klambda')
#print ("Amp R_1/2 =",pred_params_amp[1], '+/-', pred_err_amp[1], 'klambda')
#
#
#arcsec_real = convert_uv_to_lm(pred_params_real[1]*1e3)
#arcsec_err_real = (pred_err_real[1]/pred_params_real[1])*convert_uv_to_lm(pred_params_real[1]*1e3)
#print ("Real R_1/2 =", arcsec_real, '+/-', arcsec_err_real, 'arcsec')
#
#arcsec_amp = convert_uv_to_lm(pred_params_amp[1]*1e3)
#arcsec_err_amp = (pred_err_amp[1]/pred_params_amp[1])*convert_uv_to_lm(pred_params_amp[1]*1e3)
#print ("Amp R_1/2 =", arcsec_amp, '+/-', arcsec_err_amp, 'arcsec')
#
#plt.errorbar(UVDist, Real, yerr=RealError, fmt='.', color='blue', label='Real')
#plt.errorbar(UVDist, Amp, yerr=AmpError, fmt='.', color='orange', label='Amp=sqrt(Real^2+Imag^2)')
#plt.plot(np.linspace(0,100,1000), gauss1d_new(np.linspace(0,100,1000), pred_params_real[0], pred_params_real[1]), color='blue', linestyle='dashed')
#plt.plot(np.linspace(0,100,1000), gauss1d_new(np.linspace(0,100,1000), pred_params_amp[0], pred_params_amp[1]), color='orange', linestyle='dashed')
#
#
#plt.xlabel('UVdist (klambda)')
#plt.ylabel('Amp (mJy)')
##plt.ylim(bottom=-0.1)
##plt.ylim(-0.1,1)
#plt.legend()
#plt.title('Variable binwidth, first bin is sparsely populated')
##plt.savefig("fig.png")
#plt.show()
