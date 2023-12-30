# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 16:28:03 2023

@author: jansen
"""


import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.optimize import curve_fit

def ReadFITS(FileName):
    FITSFile = fits.open(FileName, lazy_load_hdu=True)
    Data = FITSFile[0].data[0][0]
    Header = FITSFile[0].header
    FITSFile.close()
    return Data, Header

def Gaussian(x, a, x0, sigma): 
    return a*np.exp(-(x-x0)**2/(2*sigma**2)) 



RMS_file = 'D:\\Master Astronomy Research year 2\\Master Project\\Profiles\\AS2COS54-NRAO-calib.split.cube.image.fits'
File = 'D:\\Master Astronomy Research year 2\\Master Project\\Profiles\\spectral_profile_AS2COS54-NRAO-calib.txt'
Cont_name = 'D:\\Master Astronomy Research year 2\\Master Project\\fig_cont_AS2COS54-NRAO-calib.png'
Profile_name = 'D:\\Master Astronomy Research year 2\\Master Project\\fig_profile_AS2COS54-NRAO-calib.png'
Title= 'AS2COS0054.1 NRAO calib'

        
Table  = open(File)
    
FITSFile = fits.open(RMS_file, lazy_load_hdu=True)
Data = FITSFile[0].data[0]
FITSFile.close()
    
RMS_cube = np.zeros(len(Data))
for i in range(len(Data)):
    RMS_cube[i] = np.nanstd(Data[i])*1e6
RMS_cube[RMS_cube<=0.0]=np.nan
RMS_cube_mean = np.nanmean(RMS_cube)
print("Mean RMS = ", round(RMS_cube_mean), "uJy/beam")

#RMS=RMS*np.sqrt(Npts/Abeam)
#RMS_mean = np.nanmean(RMS)
#print("Mean RMS = ", RMS_mean, "mJy/pix")

Vel, Flux = np.genfromtxt(Table, unpack=True)
Table.close()

Flux = Flux * 1000 #to mJy

fig = plt.figure()
ax=fig.add_subplot(111)
plt.ylabel(r'Flux density (mJy)', color = 'black', fontsize = 10)
plt.xlabel(r"Velocity offset (km/s)", fontsize = 10)

plt.tick_params(axis='both', which = 'major', length=10, direction = 'in', width=0.5, color = 'black', labelsize = 10)
plt.tick_params(axis='both', which = 'minor',length=5, direction = 'in', width=0.5, color = 'black', labelsize = 10)

y_max = np.max(Flux)
plt.ylim(-0.1,y_max+0.1)    
plt.xlim(-2000, 2000)

x_for_gaussian = Vel/1000
y_for_gaussian = Flux
x_for_model = np.linspace(-2000, 2000, 1000)/1000
      
popt, pcov = curve_fit(Gaussian, x_for_gaussian, y_for_gaussian)#, sigma=RMS_cube) 
#print(popt[0],popt[1]*100,popt[2]*100) 
  
y_for_model = Gaussian(x_for_model, popt[0], popt[1], popt[2]) 

ax.plot(x_for_model*1000, y_for_model, c='k') 

FWHM = 2*np.sqrt(2*np.log(2))*popt[2]*1000
err_FWHM = 2*np.sqrt(2*np.log(2))*np.sqrt(np.diag(pcov))[2]*1000
#print("FWHM =", round(FWHM), '+-', round(err_FWHM))

plt.step(Vel, Flux, lw = 1, where='mid', color = 'darkorange', zorder = 3)
#plt.errorbar(Vel, Flux, yerr=RMS_cube, fmt='none', capsize=4, color='gray', alpha=0.6)
#plt.fill_between(Freq, -RMS, RMS, facecolor = '0.85', edgecolor = 'none', zorder = 1)
plt.fill_between(Vel, Flux, 0, step="mid", edgecolor = 'none', zorder = 1, alpha=0.4, color = 'darkorange')

#    plt.bar(Vel-0.01, Flux, lw = 1, facecolor = 'navajowhite', edgecolor = 'none', width = 10,zorder = 2)

plt.plot([-2000,2000],[0,0], c = '0', lw = 1, linestyle = 'dashed', zorder = 3)
plt.plot([0,0],[-2,5], c = '0', lw = 1, linestyle = 'dashed', zorder = 3)

plt.text(-2000+100, y_max, Title, verticalalignment='center', fontsize = 15)
    
plt.gcf().set_size_inches(6,3)

plt.show()

#plt.savefig('spectrum_J1202_HCN32_LSB_norms.png', dpi = 200, bbox_inches='tight')
        


#%%

import numpy as np 
import matplotlib.pyplot as plt 
from scipy.optimize import curve_fit 
  
xdata = Vel/1000
ydata = Flux
      
x_for_gaussian = Vel/1000
y_for_gaussian = Flux
x_for_model = np.linspace(-2000, 2000, 1000)/1000



# Recast xdata and ydata into numpy arrays so we can use their handy features 
xdata = np.asarray(xdata) 
ydata = np.asarray(ydata) 
plt.plot(xdata, ydata, 'o') 
  
# Define the Gaussian function 
 
parameters, covariance = popt, pcov = curve_fit(Gaussian, x_for_gaussian, y_for_gaussian) 

 
  
fit_A = parameters[0] 
fit_B = parameters[1] 
fit_C = parameters[2]
  
fit_y = Gaussian(xdata, fit_A, fit_B, fit_C) 
plt.plot(xdata, ydata, 'o', label='data') 
plt.plot(xdata, fit_y, '-', label='fit') 
plt.legend()


