# -*- coding: utf-8 -*-
"""
Created on Tue Oct 10 15:30:16 2023

@author: jansen
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from matplotlib.patches import Ellipse  # this is for the beam contour
from scipy.optimize import curve_fit

#%%

def ReadFITS(FileName):
    FITSFile = fits.open(FileName, lazy_load_hdu=True)
    Data = FITSFile[0].data[0][0]
    Header = FITSFile[0].header
    FITSFile.close()
    return Data, Header

def Gaussian(x, a, x0, sigma): 
    return a*np.exp(-(x-x0)**2/(2*sigma**2)) 

def MakePlot(File, Title, Name, Save):

    Offset = np.asarray((0,0,0,0))
    Data, Header = ReadFITS(File)
    Bmaj = Header['BMAJ']*3600.0    # major axis, convert deg -> arcsec
    Bmin = Header['BMIN']*3600.0
    Bpa = Header['BPA']+90

    #print ("Beam:", Bmaj, Bmin, Bpa)

    PxScale = abs(Header['CDELT1']*3600.0) #arcsec
    
    DimX=128
    Extent0=np.asarray([-DimX*PxScale/2.,DimX*PxScale/2.,-DimX*PxScale/2.,DimX*PxScale/2.])+Offset
    
    RMS = np.std(Data[0:52,0:128])
    print ("RMS =", RMS)
    
    fig = plt.figure()
    ax=fig.add_subplot(111)
    plt.xticks([])
    plt.yticks([])
    
    plt.xlim(-10,10)
    plt.ylim(-10,10)

    ax1=plt.imshow(Data*1000., cmap = 'viridis', origin='lower', interpolation = 'none', extent =Extent0)
    plt.contour(Data/RMS, levels=[-4,-2,2,4,6,8,10,12], colors =['0'], linewidths =[1], extent =Extent0)

    cb = plt.colorbar(ax1, fraction=0.035)
    cb.set_label('mJy/beam', fontsize = 12)

    plt.text(8.25,-8.25, '2\"', color = '1',verticalalignment='center', horizontalalignment='center',backgroundcolor='none', fontsize = 16)
    plt.plot([7.25,9.25],[-9.25,-9.25], color = '1', lw = 2)
    
    Beam = Ellipse([-6.5, -6.5], width = Bmaj, height = Bmin, angle = Bpa, hatch = '/////', fc = 'none', ec = '1', lw=1, zorder = 10)
    ax.add_patch(Beam)

    plt.text(-9.25, 8.75, Title, color = '1', verticalalignment='center', fontsize = 20)
    
    if Save==True:
        fig.savefig(Name,bbox_inches='tight')
        plt.close(fig)
        
    return RMS
    
def MakeProfile(File, Titel, Name, Save):
        
    Table  = open(File)
    
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
    
    x_for_gaussian = Vel/100
    y_for_gaussian = Flux
    x_for_model = np.linspace(-2000, 2000, 1000)/100 
          
    popt, pcov = curve_fit(Gaussian, x_for_gaussian, y_for_gaussian) 
    #print(popt[0],popt[1]*100,popt[2]*100) 
      
    y_for_model = Gaussian(x_for_model, popt[0], popt[1], popt[2]) 
    
    ax.plot(x_for_model*100, y_for_model, c='k') 
    
    FWHM = 2*np.sqrt(2*np.log(2))*popt[2]*100
    print("FWHM =", FWHM)

    plt.step(Vel, Flux, lw = 1, where='mid', color = 'darkorange', zorder = 3)
    #plt.fill_between(Freq, -RMS, RMS, facecolor = '0.85', edgecolor = 'none', zorder = 1)
    plt.fill_between(Vel, Flux, 0, step="mid", edgecolor = 'none', zorder = 1, alpha=0.4, color = 'darkorange')
    
#    plt.bar(Vel-0.01, Flux, lw = 1, facecolor = 'navajowhite', edgecolor = 'none', width = 10,zorder = 2)
    
    plt.plot([-2000,2000],[0,0], c = '0', lw = 1, linestyle = 'dashed', zorder = 3)
    plt.plot([0,0],[-2,5], c = '0', lw = 1, linestyle = 'dashed', zorder = 3)
    
    plt.text(-2000+100, y_max, Title, verticalalignment='center', fontsize = 15)
        
    plt.gcf().set_size_inches(6,3)
    
    #plt.savefig('spectrum_J1202_HCN32_LSB_norms.png', dpi = 200, bbox_inches='tight')
    if Save==True:
        fig.savefig(Name,bbox_inches='tight')
        plt.close(fig)
        
    return FWHM
    

#%%
   
Savecondition = False 

#%%

Fits_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Moment_zero_images\\as2uds10_128_05_100.split.cube.image.mom0.fits'
Txt_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Profiles\\spectral_profile_as2uds10.txt'
Cont_name = 'D:\\Master Astronomy Research year 2\\Master Project\\fig_cont_AS2UDS010.0.png'
Profile_name = 'D:\\Master Astronomy Research year 2\\Master Project\\fig_profile_AS2UDS010.0.png'
Title= 'AS2UDS010.0'

RMS_AS2UD10 = MakePlot(Fits_name, Title, Cont_name, Savecondition)
FWHM_AS2UD10 = MakeProfile(Txt_name, Title, Profile_name, Savecondition)

#%%

Fits_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Moment_zero_images\\as2cos23-calib_128_05_100.split.cube.image.mom0.fits'
Txt_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Profiles\\spectral_profile_as2cos23.txt'
Cont_name = 'D:\\Master Astronomy Research year 2\\Master Project\\fig_cont_AS2COS0023.1.png'
Profile_name = 'D:\\Master Astronomy Research year 2\\Master Project\\fig_profile_AS2COS0023.1.png'
Title= 'AS2COS0023.1'

RMS_AS2COS23 = MakePlot(Fits_name, Title, Cont_name, Savecondition)
FWHM_AS2COS23 = MakeProfile(Txt_name, Title, Profile_name, Savecondition)


#%%


# perr = np.sqrt(np.diag(pcov))
# (perr[2]/popt[2])*fwhm
# Leave one out analysis Jackknife
# https://mathematica.stackexchange.com/questions/185024/error-propagation-in-finding-full-width-at-half-maximum


