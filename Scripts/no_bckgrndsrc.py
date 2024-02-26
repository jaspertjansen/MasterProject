# -*- coding: utf-8 -*-
"""
Created on Wed Feb 21 16:48:11 2024

@author: jansen
"""

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


Title = 'AS2UDS012.0'
vis = '{}_NRAO_target.ms.split.line.shifted'.format(Title)
uvbinsize = 5
UVDataFile_nobckgrndsrc = 'D:\\Master Astronomy Research year 2\\Master Project\\uvfit\\'+vis+'_binned_'+(str)(uvbinsize)+'_klambda.txt'
data_file = open(UVDataFile_nobckgrndsrc)
UVDist_no, Real_no, RealError_no, Imag_no, ImagError_no, NoVis_no = np.loadtxt(data_file, unpack = True)
data_file.close()

UVDist_no = UVDist_no[~np.isnan(UVDist_no)]
Real_no = Real_no [~np.isnan(Real_no)]
RealError_no = RealError_no [~np.isnan(RealError_no)]
Imag_no = Imag_no [~np.isnan(Imag_no)]
ImagError_no = ImagError_no [~np.isnan(ImagError_no)]

Amp_no = (Real_no**2+Imag_no**2)**(1/2)
AmpError_no = RealError_no * (Real_no/Amp_no) + ImagError_no * (Imag_no/Amp_no)




vis = '{}_NRAO_target.ms.split.line'.format(Title)
UVDataFile_bckgrndsrc = 'D:\\Master Astronomy Research year 2\\Master Project\\uvfit\\'+vis+'_binned_'+(str)(uvbinsize)+'_klambda.txt'
data_file = open(UVDataFile_bckgrndsrc)
UVDist, Real, RealError, Imag, ImagError, NoVis = np.loadtxt(data_file, unpack = True)
data_file.close()

UVDist = UVDist[~np.isnan(UVDist)]
Real = Real [~np.isnan(Real)]
RealError = RealError [~np.isnan(RealError)]
Imag = Imag [~np.isnan(Imag)]
ImagError = ImagError [~np.isnan(ImagError)]

Amp = (Real**2+Imag**2)**(1/2)
AmpError = RealError * (Real/Amp) + ImagError * (Imag/Amp)





plt.errorbar(UVDist, Real_no, yerr=RealError_no, fmt='.', label='Real noshift', color='blue')
plt.errorbar(UVDist, Real, yerr=RealError, marker='v', fmt='.', label='Real shift', color='blue')
plt.errorbar(UVDist, Amp_no, yerr=AmpError_no, label='Amp noshift', color='orange', fmt='.')
plt.errorbar(UVDist, Amp, yerr=AmpError, marker='v', label='Amp shift', color='orange', fmt='.')

plt.xlabel('UVdist (klambda)')
plt.ylabel('Amp (mJy)')
plt.legend()
plt.title('Visibility plot for UDS12 with and without phase shifting')
plt.ylim(-0.1,0.5)
plt.show()








#%%

Title = 'AS2COS0031.1'
vis = '{}_NRAO_target_nobckgrndsrc.ms.split.line.shifted'.format(Title)
uvbinsize = 5
UVDataFile_nobckgrndsrc = 'D:\\Master Astronomy Research year 2\\Master Project\\uvfit\\'+vis+'_binned_'+(str)(uvbinsize)+'_klambda.txt'
data_file = open(UVDataFile_nobckgrndsrc)
UVDist_no, Real_no, RealError_no, Imag_no, ImagError_no, NoVis_no = np.loadtxt(data_file, unpack = True)
data_file.close()

UVDist_no = UVDist_no[~np.isnan(UVDist_no)]
Real_no = Real_no [~np.isnan(Real_no)]
RealError_no = RealError_no [~np.isnan(RealError_no)]
Imag_no = Imag_no [~np.isnan(Imag_no)]
ImagError_no = ImagError_no [~np.isnan(ImagError_no)]

Amp_no = (Real_no**2+Imag_no**2)**(1/2)
AmpError_no = RealError_no * (Real_no/Amp_no) + ImagError_no * (Imag_no/Amp_no)




vis = '{}_NRAO_target.ms.split.line.shifted'.format(Title)
UVDataFile_bckgrndsrc = 'D:\\Master Astronomy Research year 2\\Master Project\\uvfit\\'+vis+'_binned_'+(str)(uvbinsize)+'_klambda.txt'
data_file = open(UVDataFile_bckgrndsrc)
UVDist, Real, RealError, Imag, ImagError, NoVis = np.loadtxt(data_file, unpack = True)
data_file.close()

UVDist = UVDist[~np.isnan(UVDist)]
Real = Real [~np.isnan(Real)]
RealError = RealError [~np.isnan(RealError)]
Imag = Imag [~np.isnan(Imag)]
ImagError = ImagError [~np.isnan(ImagError)]

Amp = (Real**2+Imag**2)**(1/2)
AmpError = RealError * (Real/Amp) + ImagError * (Imag/Amp)





plt.errorbar(UVDist, Real_no, yerr=RealError_no, fmt='.', label='Real nobckgrndsrc', color='blue')
plt.errorbar(UVDist, Real, yerr=RealError, marker='v', fmt='.', label='Real bckgrndsrc', color='blue')
plt.errorbar(UVDist, Amp_no, yerr=AmpError_no, label='Amp nobckgrndsrc', color='orange', fmt='.')
plt.errorbar(UVDist, Amp, yerr=AmpError, marker='v', label='Amp bckgrndsrc', color='orange', fmt='.')

plt.xlabel('UVdist (klambda)')
plt.ylabel('Amp (mJy)')
plt.legend()
plt.title('Visibility plot for COS31 with and without backgroundsource')
plt.ylim(-0.1,0.5)
plt.show()



#def gauss1d_new (x, a, err):
#    return a* np.exp(-(x**2)/(2*(err)**2))
#
#def straight2(x, B): # this is your 'straight line' y=f(x)
#    return 0*x+B
#
#def convert_uv_to_lm(theta_uv): #fromhttps://casadocs.readthedocs.io/en/stable/notebooks/UVTaper_Imaging_Weights.html
#    theta_lm  = 3600 * ( 4*np.log(2)/np.pi ) / (theta_uv * np.pi/180.0)
#    return theta_lm



#UVDist = UVDist[~np.isnan(UVDist)]
#Real = Real [~np.isnan(Real)]
#RealError = RealError [~np.isnan(RealError)]
#Imag = Imag [~np.isnan(Imag)]
#ImagError = ImagError [~np.isnan(ImagError)]
#
#Amp = (Real**2+Imag**2)**(1/2)
#AmpError = RealError * (Real/Amp) + ImagError * (Imag/Amp)
#
#pred_params_real, pcov_real = curve_fit(gauss1d_new, UVDist[UVDist<60], Real[UVDist<60], sigma = RealError[UVDist<60], maxfev=100000)
#pred_params_amp, pcov_amp = curve_fit(gauss1d_new, UVDist[UVDist<60], Amp[UVDist<60], sigma = AmpError[UVDist<60], maxfev=100000)
#
#pred_err_real = np.sqrt(np.diag(pcov_real))
#pred_err_amp = np.sqrt(np.diag(pcov_amp))
#
#pred_params_real2, pcov_real2 = curve_fit(straight2, UVDist[UVDist<60], Real[UVDist<60], sigma = RealError[UVDist<60], maxfev=100000)
#pred_params_amp2, pcov_amp2 = curve_fit(straight2, UVDist[UVDist<60], Amp[UVDist<60], sigma = AmpError[UVDist<60], maxfev=100000)
#
#pred_err_real2 = np.sqrt(np.diag(pcov_real2))
#pred_err_amp2 = np.sqrt(np.diag(pcov_amp2))
#
##print UVDataFile,
#print ("Real R_1/2 =",round(pred_params_real[1]), '+/-', round(pred_err_real[1]), 'klambda')
#print ("Amp R_1/2 =",round(pred_params_amp[1]), '+/-', round(pred_err_amp[1]), 'klambda')
#
#print ("Real total flux =",round(pred_params_real[0],2), '+/-', round(pred_err_real[0],2), 'mJy')
#print ("Amp total flux =",round(pred_params_amp[0],2), '+/-', round(pred_err_amp[0],2), 'mJy')
#
#print ("Real point source flux=",round(pred_params_real2[0],2), '+/-', round(pred_err_real2[0],2), 'mJy')
#print ("Amp point source flux=",round(pred_params_amp2[0],2), '+/-', round(pred_err_amp2[0],2), 'mJy')
#
#arcsec_real = convert_uv_to_lm(pred_params_real[1]*1e3)
#arcsec_err_real = (pred_err_real[1]/pred_params_real[1])*convert_uv_to_lm(pred_params_real[1]*1e3)
#print ("Real R_1/2 =", round(arcsec_real,1), '+/-', round(arcsec_err_real,1), 'arcsec')
#
#arcsec_amp = convert_uv_to_lm(pred_params_amp[1]*1e3)
#arcsec_err_amp = (pred_err_amp[1]/pred_params_amp[1])*convert_uv_to_lm(pred_params_amp[1]*1e3)
#print ("Amp R_1/2 =", round(arcsec_amp,1), '+/-', round(arcsec_err_amp,1), 'arcsec')
#
#plt.errorbar(UVDist, Real, yerr=RealError, fmt='.', color='blue')
#plt.errorbar(UVDist, Amp, yerr=AmpError, fmt='.', color='orange')
#
#plt.plot(UVDist, gauss1d_new(UVDist, pred_params_real[0], pred_params_real[1]), color='blue', linestyle='dashed', \
#         label='{} +- {} mJy, {} +- {} "'.format(round(pred_params_real[0],2), round(pred_err_real[0],2), round(arcsec_real,1),round(arcsec_err_real,1)))
#plt.plot(UVDist, gauss1d_new(UVDist, pred_params_amp[0], pred_params_amp[1]), color='orange', linestyle='dashed', \
#         label='{} +- {} mJy, {} +- {} "'.format(round(pred_params_amp[0],2), round(pred_err_amp[0],2), round(arcsec_amp,1),round(arcsec_err_amp,1)))
#
#plt.plot(UVDist, straight2(UVDist, pred_params_real2), color='blue', label='{} +- {} mJy'.format(round(pred_params_real2[0],2), round(pred_err_real2[0],2)))
#plt.plot(UVDist, straight2(UVDist, pred_params_amp2), color='orange', label='{} +- {} mJy'.format(round(pred_params_amp2[0],2), round(pred_err_amp2[0],2)))
#
#
#plt.xlabel('UVdist (klambda)')
#plt.ylabel('Amp (mJy)')
##plt.ylim(bottom=-0.1)
#plt.ylim(-0.1,0.5)
#plt.legend()
#plt.title('{}, Ignoring data after 60klambda for fit\n Blue=Real, Orange=Amp'.format(Title))
##plt.savefig('D:\\Master Astronomy Research year 2\\Master Project\\uvplanefit_{}.png'.format(Title))
#plt.show()
