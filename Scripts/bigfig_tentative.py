# -*- coding: utf-8 -*-
"""
Created on Fri Feb  2 22:28:19 2024

@author: jansen
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from matplotlib.patches import Ellipse  # this is for the beam contour
from scipy.optimize import curve_fit


def ReadFITS(FileName):
    FITSFile = fits.open(FileName, lazy_load_hdu=True)
    Data = FITSFile[0].data[0][0]
    Header = FITSFile[0].header
    FITSFile.close()
    return Data, Header

def Gaussian(x, a, x0, sigma): 
    return a*np.exp(-(x-x0)**2/(2*sigma**2)) 

def MakePlot(File, Npts_5_SNR, Abeam_SNR, Title, n):

    Offset = np.asarray((0,0,0,0))
    Data, Header = ReadFITS(File)
    Bmaj = Header['BMAJ']*3600.0    # major axis, convert deg -> arcsec
    Bmin = Header['BMIN']*3600.0
    Bpa = Header['BPA']
    
    PxScale = abs(Header['CDELT1']*3600.0) #arcsec
    
    DimX=128
    Extent0=np.asarray([-DimX*PxScale/2.,DimX*PxScale/2.,-DimX*PxScale/2.,DimX*PxScale/2.])+Offset
    
    RMS = np.std(Data[0:52,0:128])
    
    Bpa = Header['BPA'] + 90
    
    n.set_xticks([])
    n.set_yticks([])
    
    n.set_xlim(-10,10)
    n.set_ylim(-10,10)

    n.imshow(Data*1000, cmap = 'viridis', origin='lower', interpolation = 'none', extent =Extent0)
    n.contour(Data/RMS, levels=[-4,-2,2,4,6,8,10,12,14,16,18,20], colors =['0'], linewidths =[1], extent =Extent0)
    n.plot(0,0,marker="+", color="white", ms=20,zorder=10,mew=1)
    #plt.contour(Data/RMS, levels=[-4,-3,-2,2,3,4,5,6,7,8,9,10,11,12], colors =['0'], linewidths =[1], extent =Extent0)
    
    n.text(8.25,-8, '2\"', color = '1',verticalalignment='center', horizontalalignment='center',backgroundcolor='none', fontsize = 8)
    n.plot([7.25,9.25],[-9.25,-9.25], color = '1', lw = 1)

#    n.scatter(np.where(Data == np.max(Data))[0][0], np.where(Data == np.max(Data))[1][0], color='white', marker="+", s=50)
    
    
    Beam = Ellipse([-6.5, -6.5], width = Bmaj, height = Bmin, angle = Bpa, hatch = '/////////', fc = 'none', ec = '1', lw=1, zorder = 10)
    n.add_patch(Beam)

#    n.text(-9.25, 8.75, Title, color = '1', verticalalignment='center', fontsize = 20)
        
    
def MakeProfile(File, RMS_file, Npts_5_SNR, Abeam_SNR, Titel, n, Fit=True):
        
    Table  = open(File)
        
    FITSFile = fits.open(RMS_file, lazy_load_hdu=True)
    Data = FITSFile[0].data[0]
    FITSFile.close()
        
    RMS_cube = np.zeros(len(Data))
    for i in range(len(Data)):
        RMS_cube[i] = np.nanstd(Data[i])*1e6
    RMS_cube[RMS_cube<=0.0]=np.nan
    
    RMS_cube=RMS_cube*np.sqrt(Npts_5_SNR/Abeam_SNR)*1e-3
    
    Vel, Flux = np.genfromtxt(Table, unpack=True)
    Table.close()
        
    Flux = Flux * 1000 #to mJy
    y_max = np.max(Flux)
  
    if n==axs15 or n==axs35 or n==axs55 or n==axs75:
        n.set_ylabel(r'Flux (mJy)', color = 'black', fontsize = 10)
    else:
        n.yaxis.tick_right()

    if n==axs15 or n==axs25:
        n.set_xlabel(r"Velocity offset (km/s)", fontsize = 10)
        n.xaxis.set_label_position("top")
        n.xaxis.tick_top()
    else:
        n.set_xticklabels([])
        n.xaxis.tick_top()

    n.text(-2000+100, y_max, Title, verticalalignment='center', fontsize = 10, weight="bold")

    n.tick_params(axis='both', which = 'major', length=10, direction = 'in', width=0.5, color = 'black', labelsize = 10)
    n.tick_params(axis='both', which = 'minor',length=5, direction = 'in', width=0.5, color = 'black', labelsize = 10)
    
    n.set_ylim(-0.1,1.3*y_max)    
    n.set_xlim(-2000, 2000)
    
    x_for_gaussian = Vel/100
    y_for_gaussian = Flux
    x_for_model = np.linspace(-2000, 2000, 1000)/100 
          
    popt, pcov = curve_fit(Gaussian, x_for_gaussian, y_for_gaussian)
      
    y_for_model = Gaussian(x_for_model, popt[0], popt[1], popt[2]) 
    
    if Fit==True:
    
        n.plot(x_for_model*100, y_for_model, c='k') 

    n.step(Vel, Flux, lw = 1, where='mid', color = 'darkorange', zorder = 3)
    n.errorbar(Vel, Flux, yerr=RMS_cube, fmt='none', capsize=4, color='gray', alpha=0.6)
    #plt.fill_between(Freq, -RMS, RMS, facecolor = '0.85', edgecolor = 'none', zorder = 1)
    n.fill_between(Vel, Flux, 0, step="mid", edgecolor = 'none', zorder = 1, alpha=0.4, color = 'darkorange')
    
#    plt.bar(Vel-0.01, Flux, lw = 1, facecolor = 'navajowhite', edgecolor = 'none', width = 10,zorder = 2)
    
    n.plot([-2000,2000],[0,0], c = '0', lw = 1, linestyle = 'dashed', zorder = 3)
    n.plot([0,0],[-2,5], c = '0', lw = 1, linestyle = 'dashed', zorder = 3)
                        
    
fig = plt.figure(figsize=(10,5), tight_layout=True)

spec = fig.add_gridspec(4, 6)

axs1 = fig.add_subplot(spec[0, 2])
axs15 = fig.add_subplot(spec[0, 0:2])
axs2 = fig.add_subplot(spec[0, 3])
axs25 = fig.add_subplot(spec[0, 4:6])

axs3 = fig.add_subplot(spec[1, 2])
axs35 = fig.add_subplot(spec[1, 0:2])
axs4 = fig.add_subplot(spec[1, 3])
axs45 = fig.add_subplot(spec[1, 4:6])

axs5 = fig.add_subplot(spec[2, 2])
axs55 = fig.add_subplot(spec[2, 0:2])
axs6 = fig.add_subplot(spec[2, 3])
axs65 = fig.add_subplot(spec[2, 4:6])

axs7 = fig.add_subplot(spec[3, 2])
axs75 = fig.add_subplot(spec[3, 0:2])


Title= 'AS2COS0001.1'
Fits_SNR_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Moment_zero_images\\{}_NRAO_combined_target.ms.split.cube.image.SNR_mom0.fits'.format(Title)
Fits_cube_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Profiles\\{}_NRAO_combined_target.ms.split.cube.image.fits'.format(Title)
Txt_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Profiles\\spectral_profile_{}_NRAO.txt'.format(Title)
Npts_5_SNR = 81
Abeam_SNR = 68.8367
MakePlot(Fits_SNR_name, Npts_5_SNR, Abeam_SNR, Title, axs1)
MakeProfile(Txt_name, Fits_cube_name, Npts_5_SNR, Abeam_SNR, Title, axs15, Fit=False)



Title= 'AS2COS0023.1'
Fits_SNR_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Moment_zero_images\\{}_Marta_target.ms.split.contsub.cube.image.SNR_mom0.fits'.format(Title)
Fits_cube_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Profiles\\{}_Marta_target.ms.split.contsub.cube.image.fits'.format(Title)
Txt_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Profiles\\spectral_profile_{}_Marta.txt'.format(Title)
Npts_5_SNR = 81
Abeam_SNR = 63.9442
MakePlot(Fits_SNR_name, Npts_5_SNR, Abeam_SNR, Title, axs2)
MakeProfile(Txt_name, Fits_cube_name, Npts_5_SNR, Abeam_SNR, Title, axs25, Fit=False)


Title= 'AS2COS0065.1'
Fits_SNR_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Moment_zero_images\\{}_NRAO_target.ms.split.cube.image.SNR_mom0.fits'.format(Title)
Fits_cube_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Profiles\\{}_NRAO_target.ms.split.cube.image.fits'.format(Title)
Txt_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Profiles\\spectral_profile_{}_NRAO.txt'.format(Title)
Npts_5_SNR = 81
Abeam_SNR = 45.0628
MakePlot(Fits_SNR_name, Npts_5_SNR, Abeam_SNR, Title, axs3)
MakeProfile(Txt_name, Fits_cube_name, Npts_5_SNR, Abeam_SNR, Title, axs35)


Title= 'AS2UDS009.0'
Fits_SNR_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Moment_zero_images\\{}_NRAO_combined_target.ms.split.cube.image.SNR_mom0.fits'.format(Title)
Fits_cube_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Profiles\\{}_NRAO_combined_target.ms.split.cube.image.fits'.format(Title)
Txt_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Profiles\\spectral_profile_{}_NRAO.txt'.format(Title)
Npts_5_SNR = 81
Abeam_SNR = 61.1574
MakePlot(Fits_SNR_name, Npts_5_SNR, Abeam_SNR, Title, axs4)
MakeProfile(Txt_name, Fits_cube_name, Npts_5_SNR, Abeam_SNR, Title, axs45)


Title= 'AS2UDS011.0'
Fits_SNR_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Moment_zero_images\\{}_NRAO_target.ms.split.cube.image.SNR_mom0.fits'.format(Title)
Fits_cube_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Profiles\\{}_NRAO_target.ms.split.cube.image.fits'.format(Title)
Txt_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Profiles\\spectral_profile_{}_NRAO.txt'.format(Title)
Npts_5_SNR = 79
Abeam_SNR = 67.4973
MakePlot(Fits_SNR_name, Npts_5_SNR, Abeam_SNR, Title, axs5)
MakeProfile(Txt_name, Fits_cube_name, Npts_5_SNR, Abeam_SNR, Title, axs55)


Title= 'AS2UDS026.0'
Fits_SNR_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Moment_zero_images\\{}_Marta_target.ms.split.cube.image.SNR_mom0.fits'.format(Title)
Fits_cube_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Profiles\\{}_Marta_target.ms.split.cube.image.fits'.format(Title)
Txt_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Profiles\\spectral_profile_{}_Marta.txt'.format(Title)
Npts_5_SNR = 79
Abeam_SNR = 44.8977
MakePlot(Fits_SNR_name, Npts_5_SNR, Abeam_SNR, Title, axs6)
MakeProfile(Txt_name, Fits_cube_name, Npts_5_SNR, Abeam_SNR, Title, axs65, Fit=False)


Title= 'AS2UDS231.0'
Fits_SNR_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Moment_zero_images\\{}_NRAO_combined_target.ms.split.cube.image.SNR_mom0.fits'.format(Title)
Fits_cube_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Profiles\\{}_NRAO_combined_target.ms.split.cube.image.fits'.format(Title)
Txt_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Profiles\\spectral_profile_{}_NRAO.txt'.format(Title)
Npts_5_SNR = 79
Abeam_SNR = 63.4806
MakePlot(Fits_SNR_name, Npts_5_SNR, Abeam_SNR, Title, axs7)
MakeProfile(Txt_name, Fits_cube_name, Npts_5_SNR, Abeam_SNR, Title, axs75)


fig.savefig('D:\\Master Astronomy Research year 2\\Master Project\\try_bigfig_tentative.pdf')














