# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 15:31:54 2024

@author: jansen
"""

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

number_max_fit = 60 #1000

uvbinsize = 5


#full sample, median is 0.17+-0.02
list_names = ['AS2COS0008.1', 'AS2COS0013.1', 'AS2COS0014.1', 'AS2COS0028.1', 'AS2COS0031.1', 'AS2COS0054.1', 'AS2COS0066.1', \
        'AS2COS0139.1', 'AS2UDS010.0', 'AS2UDS012.0', 'AS2UDS126.0', 'AEG2', 'CDFN1', 'CDFN2', 'CDFN8']

list_vis_names = ['NRAO_combined_target.ms.split.line.shifted', 'NRAO_target.ms.split.line.shifted', 'NRAO_all_combined_target_nobckgrndsrc.ms.split.line.shifted',\
             'NRAO_target.ms.split.line.shifted', 'NRAO_target_nobckgrndsrc.ms.split.line.shifted', 'NRAO_target.ms.split.line.shifted',\
             'NRAO_combined_target.ms.split.line.shifted', 'NRAO_target.ms.split.line.shifted', 'Marta_combined_target.ms.split.line.shifted',\
             'NRAO_target.ms.split.line.shifted', 'NRAO_combined_target.ms.split.line.shifted', 'Marta_target.ms.split.line.shifted',\
             'NRAO_target.ms.split.line.shifted', 'NRAO_target.ms.split.line.shifted', 'NRAO_target_nobckgrndsrc.ms.split.line.shifted']


#selection, median is 0.19 +- 0.04
list_names = ['AS2COS0013.1', 'AS2COS0031.1', 'AS2COS0139.1', 'AS2UDS012.0', 'AEG2', 'CDFN2']

list_vis_names = ['NRAO_target.ms.split.line.shifted', 'NRAO_target_nobckgrndsrc.ms.split.line.shifted', 'NRAO_target.ms.split.line.shifted',\
                  'NRAO_target.ms.split.line.shifted', 'Marta_target.ms.split.line.shifted', 'NRAO_target.ms.split.line.shifted']


#new detections
list_names = ['AS2COS0013.1', 'AS2UDS010.0', 'AS2UDS012.0']

list_vis_names = ['NRAO_target.ms.split.2FWHM', 'Marta_combined_target.ms.split.2FWHM', 'NRAO_target.ms.split.2FWHM']

#new tentative
list_names = ['AS2COS0031.1', 'AS2COS0139.1', 'CDFN1']

list_vis_names = ['NRAO_target_nobckgrndsrc.ms.split.2FWHM', 'NRAO_target.ms.split.2FWHM', 'NRAO_target.ms.split.2FWHM']


list_UVDist = []
list_Real = []
list_RealError = []
list_Amp = []
list_AmpError = []


for i in range(len(list_names)):
    UVDataFile = 'D:\\Master Astronomy Research year 2\\Master Project\\uvfit\\'+list_names[i]+'_'+list_vis_names[i]+'_binned_'+(str)(uvbinsize)+'_klambda.txt'
    data_file = open(UVDataFile)
    UVDist, Real, RealError, Imag, ImagError, NoVis = np.loadtxt(data_file, unpack = True)
    data_file.close()
    
    UVDist = UVDist [~np.isnan(UVDist)]
    Real = Real [~np.isnan(Real)]    
    RealError = RealError [~np.isnan(RealError)]
    Imag = Imag [~np.isnan(Imag)]
    ImagError = ImagError [~np.isnan(ImagError)]
    
    Amp = (Real**2+Imag**2)**(1/2)
    AmpError = RealError * (Real/Amp) + ImagError * (Imag/Amp)
        
    missing = int(22-len(UVDist))
            
    if missing>0:
        UVDist = np.append(UVDist, np.zeros(missing) + np.nan)
        Real = np.append(Real, np.zeros(missing) + np.nan)
        RealError = np.append(RealError, np.zeros(missing) + np.nan)
        Amp = np.append(Amp, np.zeros(missing) + np.nan)
        AmpError = np.append(AmpError, np.zeros(missing) + np.nan)
        
   
    list_UVDist.append(UVDist)
    list_Real.append(Real)
    list_RealError.append(RealError)
    list_Amp.append(Amp)
    list_AmpError.append(AmpError)
    
    
mean_UVDist = np.nanmean(list_UVDist, axis=0); std_UVDist = np.nanstd(list_UVDist, axis=0)
mean_Real = np.nanmean(list_Real, axis=0); std_Real = np.nanstd(list_Real, axis=0)
#mean_RealError = np.nanmean(list_RealError, axis=0); std_RealError = np.nanstd(list_RealError, axis=0)
mean_Amp = np.nanmean(list_Amp, axis=0); std_Amp = np.nanstd(list_Amp, axis=0)
#mean_AmpError = np.nanmean(list_AmpError, axis=0); std_AmpError = np.nanstd(list_AmpError, axis=0)


#%%

        for i in range(n_it):
            random_flux_it = np.random.normal(normal_flux_median, normal_flux_err_median)
            new_flux_it = spectres(new_wave, median_wav, random_flux_it, verbose=False)
            big_array.append(new_flux_it)
            
            #Print number of iterations to keep track of the process
            if i%100 == 0:
                pass
                #print(i)
                
        big_array = np.array(big_array)
        
        #Use spectres code to define the new flux
        new_flux = spectres(new_wave, median_wav, normal_flux_median, verbose=False)
        
        #Use the Monte Carlo error calculation method to calculate the errors
        new_flux_err = np.nanstd(big_array, axis=0)










#%%


def gauss1d_new (x, a, err):
    return a* np.exp(-(x**2)/(2*(err)**2))

def straight2(x, B): # this is your 'straight line' y=f(x)
    return 0*x+B

def convert_uv_to_lm(theta_uv): #fromhttps://casadocs.readthedocs.io/en/stable/notebooks/UVTaper_Imaging_Weights.html
    theta_lm  = 3600 * ( 4*np.log(2)/np.pi ) / (theta_uv * np.pi/180.0)
    return theta_lm


pred_params_mean_Real, pcov_mean_Real = curve_fit(gauss1d_new, mean_UVDist[mean_UVDist<number_max_fit], mean_Real[mean_UVDist<number_max_fit], sigma = std_Real[mean_UVDist<number_max_fit], maxfev=100000, bounds=([0, -np.inf], [0.5, np.inf]))
pred_err_mean_Real = np.sqrt(np.diag(pcov_mean_Real))
residual_mean_Real = mean_Real - gauss1d_new(mean_UVDist, *pred_params_mean_Real)
chisq_mean_Real = sum((residual_mean_Real[mean_UVDist<number_max_fit] / std_Real[mean_UVDist<number_max_fit]) ** 2)
print("Chisq mean_Real Gaus", round(chisq_mean_Real))

pred_params_mean_Amp, pcov_mean_Amp = curve_fit(gauss1d_new, mean_UVDist[mean_UVDist<number_max_fit], mean_Amp[mean_UVDist<number_max_fit], sigma = std_Amp[mean_UVDist<number_max_fit], maxfev=100000, bounds=([0, -np.inf], [0.5, np.inf]))
pred_err_mean_Amp = np.sqrt(np.diag(pcov_mean_Amp))
residual_mean_Amp = mean_Amp - gauss1d_new(mean_UVDist, *pred_params_mean_Amp)
chisq_mean_Amp = sum((residual_mean_Amp[mean_UVDist<number_max_fit] / std_Amp[mean_UVDist<number_max_fit]) ** 2)
print("Chisq mean_Amp Gaus", round(chisq_mean_Amp))

pred_params_mean_Real2, pcov_mean_Real2 = curve_fit(straight2, mean_UVDist[mean_UVDist<number_max_fit], mean_Real[mean_UVDist<number_max_fit], sigma = std_Real[mean_UVDist<number_max_fit], maxfev=100000, bounds=([0], [0.5]))
pred_err_mean_Real2 = np.sqrt(np.diag(pcov_mean_Real2))
residual_mean_Real2 = mean_Real - straight2(mean_UVDist, *pred_params_mean_Real2)
chisq_mean_Real2 = sum((residual_mean_Real2[mean_UVDist<number_max_fit] / std_Real[mean_UVDist<number_max_fit]) ** 2)
print("Chisq mean_Real point source", round(chisq_mean_Real2))

pred_params_mean_Amp2, pcov_mean_Amp2 = curve_fit(straight2, mean_UVDist[mean_UVDist<number_max_fit], mean_Amp[mean_UVDist<number_max_fit], sigma = std_Amp[mean_UVDist<number_max_fit], maxfev=100000, bounds=([0], [0.5]))
pred_err_mean_Amp2 = np.sqrt(np.diag(pcov_mean_Amp2))
residual_mean_Amp2 = mean_Amp - straight2(mean_UVDist, *pred_params_mean_Amp2)
chisq_mean_Amp2 = sum((residual_mean_Amp2[mean_UVDist<number_max_fit] / std_Amp[mean_UVDist<number_max_fit]) ** 2)
print("Chisq mean_Amp point source", round(chisq_mean_Amp2))


print ("mean_Real R_1/2 =",round(pred_params_mean_Real[1]), '+/-', round(pred_err_mean_Real[1]), 'klambda')
print ("mean_Amp R_1/2 =",round(pred_params_mean_Amp[1]), '+/-', round(pred_err_mean_Amp[1]), 'klambda')

arcsec_mean_Real = convert_uv_to_lm(pred_params_mean_Real[1]*1e3)
arcsec_err_mean_Real = (pred_err_mean_Real[1]/pred_params_mean_Real[1])*convert_uv_to_lm(pred_params_mean_Real[1]*1e3)
print ("mean_Real R_1/2 =", round(arcsec_mean_Real,1), '+/-', round(arcsec_err_mean_Real,1), 'arcsec')

arcsec_mean_Amp = convert_uv_to_lm(pred_params_mean_Amp[1]*1e3)
arcsec_err_mean_Amp = (pred_err_mean_Amp[1]/pred_params_mean_Amp[1])*convert_uv_to_lm(pred_params_mean_Amp[1]*1e3)
print ("mean_Amp R_1/2 =", round(arcsec_mean_Amp,1), '+/-', round(arcsec_err_mean_Amp,1), 'arcsec')


print ("mean_Real Gaus total flux =",round(pred_params_mean_Real[0],2), '+/-', round(pred_err_mean_Real[0],2), 'mJy')
print ("mean_Amp Gaus total flux =",round(pred_params_mean_Amp[0],2), '+/-', round(pred_err_mean_Amp[0],2), 'mJy')

print ("mean_Real point source total flux=",round(pred_params_mean_Real2[0],2), '+/-', round(pred_err_mean_Real2[0],2), 'mJy')
print ("mean_Amp point source total flux=",round(pred_params_mean_Amp2[0],2), '+/-', round(pred_err_mean_Amp2[0],2), 'mJy')


plot_UVDist = np.linspace(min(mean_UVDist), max(mean_UVDist),10000)


plt.plot(plot_UVDist, gauss1d_new(plot_UVDist, pred_params_mean_Real[0], pred_params_mean_Real[1]), color='blue', linestyle='dashed', 
         label='chisq={}, {} +- {} mJy, {} +- {} "'.format(round(chisq_mean_Real), round(pred_params_mean_Real[0],2), round(pred_err_mean_Real[0],2), round(arcsec_mean_Real,1),round(arcsec_err_mean_Real,1)))
plt.plot(plot_UVDist, gauss1d_new(plot_UVDist, pred_params_mean_Amp[0], pred_params_mean_Amp[1]), color='orange', linestyle='dashed', 
         label='chisq={}, {} +- {} mJy, {} +- {} "'.format(round(chisq_mean_Amp), round(pred_params_mean_Amp[0],2), round(pred_err_mean_Amp[0],2), round(arcsec_mean_Amp,1),round(arcsec_err_mean_Amp,1)))

plt.plot(plot_UVDist, straight2(plot_UVDist, pred_params_mean_Real2), color='blue', label='chisq={}, {} +- {} mJy'.format(round(chisq_mean_Real2), round(pred_params_mean_Real2[0],2), round(pred_err_mean_Real2[0],2)))
plt.plot(plot_UVDist, straight2(plot_UVDist, pred_params_mean_Amp2), color='orange', label='chisq={}, {} +- {} mJy'.format(round(chisq_mean_Amp2), round(pred_params_mean_Amp2[0],2), round(pred_err_mean_Amp2[0],2)))

plt.plot([0,max(UVDist)],[0,0], c = '0', lw = 1, linestyle = 'dashed', zorder = 3)


plt.xlabel('UVdist (klambda)')
plt.ylabel('Amp (mJy)')
#plt.ylim(bottom=-0.1)
plt.ylim(-0.25,0.75)
plt.xlim(0,max(UVDist))
plt.legend()

plt.errorbar(mean_UVDist, mean_Real, yerr=std_Real, fmt='.', color='blue')
plt.errorbar(mean_UVDist, mean_Amp, yerr=std_Amp, fmt='.', color='orange')
plt.title('Poor mans stack, no weighting used\n Blue=Real, Orange=Amp')
#plt.savefig('D:\\Master Astronomy Research year 2\\Master Project\\poormansstack_nobounds.png')


plt.show()
    

