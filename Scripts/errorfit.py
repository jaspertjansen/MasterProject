import numpy as np
import math
import scipy as sp
from scipy.optimize import curve_fit

import matplotlib.pyplot as plt
import sys

# Path to the file with binned uv data
#UVDataFile =     # e.g. data_binned.txt



UVDataFile = vis+'_binned_'+(str)(uvbinsize)+'_klambda.txt'

def gauss1d_new (x, a, err):
    return a* np.exp(-(x**2)/(2*(err)**2))

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

pred_params_real, pcov_real = curve_fit(gauss1d_new, UVDist[UVDist<60], Real[UVDist<60], sigma = RealError[UVDist<60], p0= [0.15,30], maxfev=100000)
pred_params_amp, pcov_amp = curve_fit(gauss1d_new, UVDist[UVDist<60], Amp[UVDist<60], sigma = AmpError[UVDist<60], p0= [0.15,30], maxfev=100000)

pred_err_real = np.sqrt(np.diag(pcov_real))
pred_err_amp = np.sqrt(np.diag(pcov_amp))


#print UVDataFile,
print ("Real R_1/2 =",pred_params_real[1], '+/-', pred_err_real[1], 'klambda')
print ("Amp R_1/2 =",pred_params_amp[1], '+/-', pred_err_amp[1], 'klambda')

arcsec_real = convert_uv_to_lm(pred_params_real[1]*1e3)
arcsec_err_real = (pred_err_real[1]/pred_params_real[1])*convert_uv_to_lm(pred_params_real[1]*1e3)
print ("Real R_1/2 =", arcsec_real, '+/-', arcsec_err_real, 'arcsec')

arcsec_amp = convert_uv_to_lm(pred_params_amp[1]*1e3)
arcsec_err_amp = (pred_err_amp[1]/pred_params_amp[1])*convert_uv_to_lm(pred_params_amp[1]*1e3)
print ("Amp R_1/2 =", arcsec_amp, '+/-', arcsec_err_amp, 'arcsec')

plt.errorbar(UVDist, Real, yerr=RealError, fmt='.', color='blue', label='Real')
plt.errorbar(UVDist, Amp, yerr=AmpError, fmt='.', color='orange', label='Amp=sqrt(Real^2+Imag^2)')
plt.plot(UVDist, gauss1d_new(UVDist, pred_params_real[0], pred_params_real[1]), color='blue', linestyle='dashed')
plt.plot(UVDist, gauss1d_new(UVDist, pred_params_amp[0], pred_params_amp[1]), color='orange', linestyle='dashed')

plt.xlabel('UVdist (klambda)')
plt.ylabel('Amp (mJy)')
#plt.ylim(bottom=-0.1)
#plt.ylim(-0.1,0.5)
plt.legend()
plt.title('Ignoring data after 60klambda for fit, if included fit for Amp does not work, \n but for Real still works')
#plt.savefig("fig.png")
plt.show()