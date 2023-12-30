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

def MakePlot(File, Npts_5_SNR, Abeam_SNR, Title, Name, Save):

    Offset = np.asarray((0,0,0,0))
    Data, Header = ReadFITS(File)
    Bmaj = Header['BMAJ']*3600.0    # major axis, convert deg -> arcsec
    Bmin = Header['BMIN']*3600.0
    Bpa = Header['BPA']
    
    #print ("Beam:", Bmaj, Bmin, Bpa)

    PxScale = abs(Header['CDELT1']*3600.0) #arcsec
    
    DimX=128
    Extent0=np.asarray([-DimX*PxScale/2.,DimX*PxScale/2.,-DimX*PxScale/2.,DimX*PxScale/2.])+Offset
    
    RMS = np.std(Data[0:52,0:128])
    #RMS = RMS * np.sqrt(Npts_5_SNR/Abeam_SNR)
    #print("RMS =", RMS)
    print("maj x min, pa in arcsec and degrees", round(Bmaj,1), round(Bmin,1), round(Bpa))
    
    Bpa = Header['BPA'] + 90
    
    fig = plt.figure()
    ax=fig.add_subplot(111)
    plt.xticks([])
    plt.yticks([])
    
    plt.xlim(-10,10)
    plt.ylim(-10,10)

    ax1=plt.imshow(Data*1000, cmap = 'viridis', origin='lower', interpolation = 'none', extent =Extent0)
    plt.contour(Data/RMS, levels=[-4,-2,2,4,6,8,10,12], colors =['0'], linewidths =[1], extent =Extent0)
    #plt.contour(Data/RMS, levels=[-4,-3,-2,2,3,4,5,6,7,8,9,10,11,12], colors =['0'], linewidths =[1], extent =Extent0)
    
    #print(RMS, Data[64,64], Data[64,64]/RMS)

    cb = plt.colorbar(ax1, fraction=0.035)
    cb.set_label('mJy/beam', fontsize = 12)

    plt.text(8.25,-8.25, '2\"', color = '1',verticalalignment='center', horizontalalignment='center',backgroundcolor='none', fontsize = 16)
    plt.plot([7.25,9.25],[-9.25,-9.25], color = '1', lw = 2)
    plt.scatter(np.where(Data == np.max(Data))[0][0], np.where(Data == np.max(Data))[1][0], color='white', marker="+", s=50)
    
#    print(Extent0)
#    
#    i,j = np.unravel_index(Data.argmax(), Data.shape)
#    print(i,j)
#    plt.scatter(j-64.25, i-63.75, color='red', marker="+", s=100)
    
    Beam = Ellipse([-6.5, -6.5], width = Bmaj, height = Bmin, angle = Bpa, hatch = '/////', fc = 'none', ec = '1', lw=1, zorder = 10)
    ax.add_patch(Beam)

    plt.text(-9.25, 8.75, Title, color = '1', verticalalignment='center', fontsize = 20)
    
    if Save==True:
        fig.savefig(Name,bbox_inches='tight')
        plt.close(fig)
        
    
def MakeProfile(File, RMS_file, Npts_5_SNR, Abeam_SNR, Titel, Name, Save):
        
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
    
    RMS_cube=RMS_cube*np.sqrt(Npts_5_SNR/Abeam_SNR)*1e-3
    RMS_cube_mean = np.nanmean(RMS_cube)
    
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
    err_FWHM = 2*np.sqrt(2*np.log(2))*np.sqrt(np.diag(pcov))[2]*100
    print("FWHM =", round(FWHM), '+-', round(err_FWHM))

    plt.step(Vel, Flux, lw = 1, where='mid', color = 'darkorange', zorder = 3)
    plt.errorbar(Vel, Flux, yerr=RMS_cube, fmt='none', capsize=4, color='gray', alpha=0.6)
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
        
    
def calc_param(Fits_SNR, Fits_Ico, Npts_5_SNR, Npts_5_Ico, Npts_12_Ico, SNR_5, Ico_5, Ico_12, Abeam_SNR, Abeam_Ico, z, D_L):
    
    Data, Header = ReadFITS(Fits_Ico)
    RMS_Ico = np.std(Data[0:52,0:128])

    
    Data, Header = ReadFITS(Fits_SNR)
    RMS_SNR = np.std(Data[0:52,0:128])
    
    RMS_SNR = RMS_SNR * np.sqrt(Npts_5_SNR/Abeam_SNR)
    print("SNR =", round(SNR_5/RMS_SNR,1))
    
        
    if SNR_5/RMS_SNR > 3:
        RMS_Ico_12 = RMS_Ico * np.sqrt(Npts_12_Ico/Abeam_Ico)
        print("I_co = ", round(Ico_12,2), "+-", round(RMS_Ico_12,2))
        L_co = calc_L(Ico_12,z, D_L)      
        print("L_co = ", round(L_co*10**(-10),1), "+-", round(L_co*10**(-10)*(RMS_Ico_12/Ico_12),1))

    else:
        RMS_Ico_5 = RMS_Ico * np.sqrt(Npts_5_Ico/Abeam_Ico)
        print("I_co = ", round(1.8*Ico_5,2), "+-", round(1.8*RMS_Ico_5,2))
        L_co = calc_L(1.8*Ico_5,z, D_L)      
        print("L_co = ", round(L_co*10**(-10),1), "+-", round(L_co*10**(-10)*(RMS_Ico_5/Ico_5),1))

      

def calc_L(Ico_12,z,D_L):
    
    v_obs = 115.271203/(1+z)
    
    return 3.25e7 * Ico_12 * (1+z)**(-3) * v_obs**(-2) * D_L**2

#%%
   
Savecondition = True

#%%

Title= 'AS2COS0014.1'
Cont_name = 'D:\\Master Astronomy Research year 2\\Master Project\\fig_cont_{}_NRAO.png'.format(Title)
Profile_name = 'D:\\Master Astronomy Research year 2\\Master Project\\fig_profile_{}_NRAO.png'.format(Title)

Fits_SNR_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Moment_zero_images\\{}_NRAO_combined_target.ms.split.cube.image.SNR_mom0.fits'.format(Title)
Fits_Ico_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Moment_zero_images\\{}_NRAO_combined_target.ms.split.cube.image.Ico_mom0.fits'.format(Title)
Fits_cube_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Profiles\\{}_NRAO_combined_target.ms.split.cube.image.fits'.format(Title)
Txt_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Profiles\\spectral_profile_{}_NRAO.txt'.format(Title)

z = 2.921
D_L = 25168.5       # Mpc, use https://www.astro.ucla.edu/%7Ewright/CosmoCalc.html with H0 = 67.8, Om = 0.31 and Ovac = 0.69

Npts_5_SNR = 81
Npts_5_Ico = 81
Npts_12_Ico = 441
SNR_5 = 1.3505e-1
Ico_5 = 1.0465e-1
Ico_12 = 2.9032e-1
Abeam_SNR = 39.4964
Abeam_Ico = 39.8679


MakePlot(Fits_SNR_name, Npts_5_SNR, Abeam_SNR, Title, Cont_name, Savecondition)
MakeProfile(Txt_name, Fits_cube_name, Npts_5_SNR, Abeam_SNR, Title, Profile_name, Savecondition)
calc_param(Fits_SNR_name, Fits_Ico_name, Npts_5_SNR, Npts_5_Ico, Npts_12_Ico, SNR_5, Ico_5, Ico_12, Abeam_SNR, Abeam_Ico, z, D_L)






#%%

Title= 'AS2UDS072.0'
Cont_name = 'D:\\Master Astronomy Research year 2\\Master Project\\fig_cont_{}_NRAO.png'.format(Title)
Profile_name = 'D:\\Master Astronomy Research year 2\\Master Project\\fig_profile_{}_NRAO.png'.format(Title)

Fits_SNR_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Moment_zero_images\\{}_NRAO_combined_no_statwt_target.ms.split.contsub.cube.image.SNR_mom0.fits'.format(Title)
Fits_Ico_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Moment_zero_images\\{}_NRAO_combined_no_statwt_target.ms.split.contsub.cube.image.Ico_mom0.fits'.format(Title)
Fits_cube_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Profiles\\{}_NRAO_combined_no_statwt_target.ms.split.contsub.cube.image.fits'.format(Title)
Txt_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Profiles\\spectral_profile_{}_NRAO.txt'.format(Title)

z = 2.406
D_L = 19915.3       # Mpc, use https://www.astro.ucla.edu/%7Ewright/CosmoCalc.html with H0 = 67.8, Om = 0.31 and Ovac = 0.69

Npts_5_SNR = 79
Npts_5_Ico = 79
Npts_12_Ico = 439
SNR_5 = 6.4494e-3
Ico_5 = -1.7931e-2
Ico_12 = 8.0721e-2
Abeam_SNR = 33.525
Abeam_Ico = 33.5464


MakePlot(Fits_SNR_name, Npts_5_SNR, Abeam_SNR, Title, Cont_name, Savecondition)
MakeProfile(Txt_name, Fits_cube_name, Npts_5_SNR, Abeam_SNR, Title, Profile_name, Savecondition)
calc_param(Fits_SNR_name, Fits_Ico_name, Npts_5_SNR, Npts_5_Ico, Npts_12_Ico, SNR_5, Ico_5, Ico_12, Abeam_SNR, Abeam_Ico, z, D_L)




#%%


Title= 'AS2UDS009.0'
Cont_name = 'D:\\Master Astronomy Research year 2\\Master Project\\fig_cont_{}_NRAO.png'.format(Title)
Profile_name = 'D:\\Master Astronomy Research year 2\\Master Project\\fig_profile_{}_NRAO.png'.format(Title)

Fits_SNR_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Moment_zero_images\\{}_NRAO_combined_target.ms.split.cube.image.SNR_mom0.fits'.format(Title)
Fits_Ico_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Moment_zero_images\\{}_NRAO_combined_target.ms.split.cube.image.Ico_mom0.fits'.format(Title)
Fits_cube_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Profiles\\{}_NRAO_combined_target.ms.split.cube.image.fits'.format(Title)
Txt_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Profiles\\spectral_profile_{}_NRAO.txt'.format(Title)

z = 2.942
D_L = 25386.0       # Mpc, use https://www.astro.ucla.edu/%7Ewright/CosmoCalc.html with H0 = 67.8, Om = 0.31 and Ovac = 0.69

Npts_5_SNR = 81
Npts_5_Ico = 79
Npts_12_Ico = 437
SNR_5 = 5.7561e-2
Ico_5 = 6.8153e-2
Ico_12 = 2.3201e-1
Abeam_SNR = 61.1574
Abeam_Ico = 61.2763


MakePlot(Fits_SNR_name, Npts_5_SNR, Abeam_SNR, Title, Cont_name, Savecondition)
MakeProfile(Txt_name, Fits_cube_name, Npts_5_SNR, Abeam_SNR, Title, Profile_name, Savecondition)
calc_param(Fits_SNR_name, Fits_Ico_name, Npts_5_SNR, Npts_5_Ico, Npts_12_Ico, SNR_5, Ico_5, Ico_12, Abeam_SNR, Abeam_Ico, z, D_L)




#%%


Title= 'AS2UDS231.0'
Cont_name = 'D:\\Master Astronomy Research year 2\\Master Project\\fig_cont_{}_NRAO.png'.format(Title)
Profile_name = 'D:\\Master Astronomy Research year 2\\Master Project\\fig_profile_{}_NRAO.png'.format(Title)

Fits_SNR_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Moment_zero_images\\{}_NRAO_combined_target.ms.split.contsub.cube.image.SNR_mom0.fits'.format(Title)
Fits_Ico_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Moment_zero_images\\{}_NRAO_combined_target.ms.split.contsub.cube.image.Ico_mom0.fits'.format(Title)
Fits_cube_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Profiles\\{}_NRAO_combined_target.ms.split.contsub.cube.image.fits'.format(Title)
Txt_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Profiles\\spectral_profile_{}_NRAO.txt'.format(Title)

z = 3.119
D_L = 27229.4       # Mpc, use https://www.astro.ucla.edu/%7Ewright/CosmoCalc.html with H0 = 67.8, Om = 0.31 and Ovac = 0.69

Npts_5_SNR = 79
Npts_5_Ico = 79
Npts_12_Ico = 439
SNR_5 = 8.8303e-2
Ico_5 = 1.0833e-1
Ico_12 = 2.0925e-1
Abeam_SNR = 63.4806
Abeam_Ico = 63.5617


MakePlot(Fits_SNR_name, Npts_5_SNR, Abeam_SNR, Title, Cont_name, Savecondition)
MakeProfile(Txt_name, Fits_cube_name, Npts_5_SNR, Abeam_SNR, Title, Profile_name, Savecondition)
calc_param(Fits_SNR_name, Fits_Ico_name, Npts_5_SNR, Npts_5_Ico, Npts_12_Ico, SNR_5, Ico_5, Ico_12, Abeam_SNR, Abeam_Ico, z, D_L)





#%%

Title= 'AS2COS0044.1'
Cont_name = 'D:\\Master Astronomy Research year 2\\Master Project\\fig_cont_{}_NRAO.png'.format(Title)
Profile_name = 'D:\\Master Astronomy Research year 2\\Master Project\\fig_profile_{}_NRAO.png'.format(Title)

Fits_SNR_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Moment_zero_images\\{}_NRAO_target.ms.split.cube.image.SNR_mom0.fits'.format(Title)
Fits_Ico_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Moment_zero_images\\{}_NRAO_target.ms.split.cube.image.Ico_mom0.fits'.format(Title)
Fits_cube_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Profiles\\{}_NRAO_target.ms.split.cube.image.fits'.format(Title)
Txt_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Profiles\\spectral_profile_{}_NRAO.txt'.format(Title)

z = 2.580
D_L = 21671.3       # Mpc, use https://www.astro.ucla.edu/%7Ewright/CosmoCalc.html with H0 = 67.8, Om = 0.31 and Ovac = 0.69

Npts_5_SNR = 81
Npts_5_Ico = 81
Npts_12_Ico = 441
SNR_5 = 1.0120e-2
Ico_5 = 4.4606e-2
Ico_12 = 2.0597e-2
Abeam_SNR = 41.5288
Abeam_Ico = 41.7937


MakePlot(Fits_SNR_name, Npts_5_SNR, Abeam_SNR, Title, Cont_name, Savecondition)
MakeProfile(Txt_name, Fits_cube_name, Npts_5_SNR, Abeam_SNR, Title, Profile_name, Savecondition)
calc_param(Fits_SNR_name, Fits_Ico_name, Npts_5_SNR, Npts_5_Ico, Npts_12_Ico, SNR_5, Ico_5, Ico_12, Abeam_SNR, Abeam_Ico, z, D_L)


#%%



Title= 'AS2COS0028.1'
Cont_name = 'D:\\Master Astronomy Research year 2\\Master Project\\fig_cont_{}_NRAO.png'.format(Title)
Profile_name = 'D:\\Master Astronomy Research year 2\\Master Project\\fig_profile_{}_NRAO.png'.format(Title)

Fits_SNR_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Moment_zero_images\\{}_NRAO_target.ms.split.cube.image.SNR_mom0.fits'.format(Title)
Fits_Ico_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Moment_zero_images\\{}_NRAO_target.ms.split.cube.image.Ico_mom0.fits'.format(Title)
Fits_cube_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Profiles\\{}_NRAO_target.ms.split.cube.image.fits'.format(Title)
Txt_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Profiles\\spectral_profile_{}_NRAO.txt'.format(Title)

z = 3.097
D_L = 26999.4       # Mpc, use https://www.astro.ucla.edu/%7Ewright/CosmoCalc.html with H0 = 67.8, Om = 0.31 and Ovac = 0.69

Npts_5_SNR = 81
Npts_5_Ico = 81
Npts_12_Ico = 441
SNR_5 = 1.3063e-1
Ico_5 = 1.5696e-1
Ico_12 = 1.5100e-1
Abeam_SNR = 41.7503
Abeam_Ico = 41.8256


MakePlot(Fits_SNR_name, Npts_5_SNR, Abeam_SNR, Title, Cont_name, Savecondition)
MakeProfile(Txt_name, Fits_cube_name, Npts_5_SNR, Abeam_SNR, Title, Profile_name, Savecondition)
calc_param(Fits_SNR_name, Fits_Ico_name, Npts_5_SNR, Npts_5_Ico, Npts_12_Ico, SNR_5, Ico_5, Ico_12, Abeam_SNR, Abeam_Ico, z, D_L)


#%%

Title= 'AS2COS0008.1'
Cont_name = 'D:\\Master Astronomy Research year 2\\Master Project\\fig_cont_{}_NRAO.png'.format(Title)
Profile_name = 'D:\\Master Astronomy Research year 2\\Master Project\\fig_profile_{}_NRAO.png'.format(Title)

Fits_SNR_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Moment_zero_images\\{}_NRAO_combined_target.ms.split.contsub.cube.image.SNR_mom0.fits'.format(Title)
Fits_Ico_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Moment_zero_images\\{}_NRAO_combined_target.ms.split.contsub.cube.image.Ico_mom0.fits'.format(Title)
Fits_cube_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Profiles\\{}_NRAO_combined_target.ms.split.contsub.cube.image.fits'.format(Title)
Txt_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Profiles\\spectral_profile_{}_NRAO.txt'.format(Title)

z = 3.581
D_L = 32113.4       # Mpc, use https://www.astro.ucla.edu/%7Ewright/CosmoCalc.html with H0 = 67.8, Om = 0.31 and Ovac = 0.69

Npts_5_SNR = 81
Npts_5_Ico = 81
Npts_12_Ico = 441
SNR_5 = 8.428e-2
Ico_5 = 1.0104e-1
Ico_12 = 1.3761e-1
Abeam_SNR = 56.2641
Abeam_Ico = 56.3711


MakePlot(Fits_SNR_name, Npts_5_SNR, Abeam_SNR, Title, Cont_name, Savecondition)
MakeProfile(Txt_name, Fits_cube_name, Npts_5_SNR, Abeam_SNR, Title, Profile_name, Savecondition)
calc_param(Fits_SNR_name, Fits_Ico_name, Npts_5_SNR, Npts_5_Ico, Npts_12_Ico, SNR_5, Ico_5, Ico_12, Abeam_SNR, Abeam_Ico, z, D_L)



#%%


Title= 'AS2UDS011.0'
Cont_name = 'D:\\Master Astronomy Research year 2\\Master Project\\fig_cont_{}_NRAO.png'.format(Title)
Profile_name = 'D:\\Master Astronomy Research year 2\\Master Project\\fig_profile_{}_NRAO.png'.format(Title)

Fits_SNR_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Moment_zero_images\\{}_NRAO_target.ms.split.cube.image.SNR_mom0.fits'.format(Title)
Fits_Ico_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Moment_zero_images\\{}_NRAO_target.ms.split.cube.image.Ico_mom0.fits'.format(Title)
Fits_cube_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Profiles\\{}_NRAO_target.ms.split.cube.image.fits'.format(Title)
Txt_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Profiles\\spectral_profile_{}_NRAO.txt'.format(Title)

z = 4.073
D_L = 37414.8       # Mpc, use https://www.astro.ucla.edu/%7Ewright/CosmoCalc.html with H0 = 67.8, Om = 0.31 and Ovac = 0.69

Npts_5_SNR = 79
Npts_5_Ico = 79
Npts_12_Ico = 439
SNR_5 = 4.1136e-2
Ico_5 = 5.8689e-2
Ico_12 = 1.3132e-1
Abeam_SNR = 67.4973
Abeam_Ico = 68.2418


MakePlot(Fits_SNR_name, Npts_5_SNR, Abeam_SNR, Title, Cont_name, Savecondition)
MakeProfile(Txt_name, Fits_cube_name, Npts_5_SNR, Abeam_SNR, Title, Profile_name, Savecondition)
calc_param(Fits_SNR_name, Fits_Ico_name, Npts_5_SNR, Npts_5_Ico, Npts_12_Ico, SNR_5, Ico_5, Ico_12, Abeam_SNR, Abeam_Ico, z, D_L)




#%%


Title= 'AS2UDS010.0'
Cont_name = 'D:\\Master Astronomy Research year 2\\Master Project\\fig_cont_{}_Marta.png'.format(Title)
Profile_name = 'D:\\Master Astronomy Research year 2\\Master Project\\fig_profile_{}_Marta.png'.format(Title)

Fits_SNR_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Moment_zero_images\\{}_Marta_combined_target.ms.split.cube.image.SNR_mom0.fits'.format(Title)
Fits_Ico_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Moment_zero_images\\{}_Marta_combined_target.ms.split.cube.image.Ico_mom0.fits'.format(Title)
Fits_cube_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Profiles\\{}_Marta_combined_target.ms.split.cube.image.fits'.format(Title)
Txt_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Profiles\\spectral_profile_{}_Marta.txt'.format(Title)

z = 3.169
D_L = 27753.0       # Mpc, use https://www.astro.ucla.edu/%7Ewright/CosmoCalc.html with H0 = 67.8, Om = 0.31 and Ovac = 0.69

Npts_5_SNR = 79
Npts_5_Ico = 79
Npts_12_Ico = 439
SNR_5 = 1.2578e-1
Ico_5 = 1.3622e-1
Ico_12 = 3.4213e-1
Abeam_SNR = 54.9011
Abeam_Ico = 54.9732


MakePlot(Fits_SNR_name, Npts_5_SNR, Abeam_SNR, Title, Cont_name, Savecondition)
MakeProfile(Txt_name, Fits_cube_name, Npts_5_SNR, Abeam_SNR, Title, Profile_name, Savecondition)
calc_param(Fits_SNR_name, Fits_Ico_name, Npts_5_SNR, Npts_5_Ico, Npts_12_Ico, SNR_5, Ico_5, Ico_12, Abeam_SNR, Abeam_Ico, z, D_L)



#%%


Title= 'AS2COS0011.1'
Cont_name = 'D:\\Master Astronomy Research year 2\\Master Project\\fig_cont_{}_NRAO.png'.format(Title)
Profile_name = 'D:\\Master Astronomy Research year 2\\Master Project\\fig_profile_{}_NRAO.png'.format(Title)

Fits_SNR_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Moment_zero_images\\{}_NRAO_target.ms.split.contsub.cube.image.SNR_mom0.fits'.format(Title)
Fits_Ico_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Moment_zero_images\\{}_NRAO_target.ms.split.contsub.cube.image.Ico_mom0.fits'.format(Title)
Fits_cube_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Profiles\\{}_NRAO_target.ms.split.contsub.cube.image.fits'.format(Title)
Txt_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Profiles\\spectral_profile_{}_NRAO.txt'.format(Title)

z = 4.786
D_L = 45249.9       # Mpc, use https://www.astro.ucla.edu/%7Ewright/CosmoCalc.html with H0 = 67.8, Om = 0.31 and Ovac = 0.69

Npts_5_SNR = 81
Npts_5_Ico = 81
Npts_12_Ico = 441
SNR_5 = 2.0673e-2
Ico_5 = 4.5381e-2
Ico_12 = -4.9697e-2
Abeam_SNR = 46.6536
Abeam_Ico = 48.3674


MakePlot(Fits_SNR_name, Npts_5_SNR, Abeam_SNR, Title, Cont_name, Savecondition)
MakeProfile(Txt_name, Fits_cube_name, Npts_5_SNR, Abeam_SNR, Title, Profile_name, Savecondition)
calc_param(Fits_SNR_name, Fits_Ico_name, Npts_5_SNR, Npts_5_Ico, Npts_12_Ico, SNR_5, Ico_5, Ico_12, Abeam_SNR, Abeam_Ico, z, D_L)


#%%

Title= 'AS2UDS012.0'
Cont_name = 'D:\\Master Astronomy Research year 2\\Master Project\\fig_cont_{}_NRAO.png'.format(Title)
Profile_name = 'D:\\Master Astronomy Research year 2\\Master Project\\fig_profile_{}_NRAO.png'.format(Title)

Fits_SNR_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Moment_zero_images\\{}_NRAO_target.ms.split.cube.image.SNR_mom0.fits'.format(Title)
Fits_Ico_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Moment_zero_images\\{}_NRAO_target.ms.split.cube.image.Ico_mom0.fits'.format(Title)
Fits_cube_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Profiles\\{}_NRAO_target.ms.split.cube.image.fits'.format(Title)
Txt_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Profiles\\spectral_profile_{}_NRAO.txt'.format(Title)

z = 2.52
D_L = 21063.5       # Mpc, use https://www.astro.ucla.edu/%7Ewright/CosmoCalc.html with H0 = 67.8, Om = 0.31 and Ovac = 0.69

Npts_5_SNR = 79
Npts_5_Ico = 79
Npts_12_Ico = 439
SNR_5 = 9.3226e-2
Ico_5 = 1.4552e-1
Ico_12 = 2.1825e-1
Abeam_SNR = 35.0724
Abeam_Ico = 35.3822


MakePlot(Fits_SNR_name, Npts_5_SNR, Abeam_SNR, Title, Cont_name, Savecondition)
MakeProfile(Txt_name, Fits_cube_name, Npts_5_SNR, Abeam_SNR, Title, Profile_name, Savecondition)
calc_param(Fits_SNR_name, Fits_Ico_name, Npts_5_SNR, Npts_5_Ico, Npts_12_Ico, SNR_5, Ico_5, Ico_12, Abeam_SNR, Abeam_Ico, z, D_L)



#%%


Title= 'AS2COS0023.1'
Cont_name = 'D:\\Master Astronomy Research year 2\\Master Project\\fig_cont_{}_Marta.png'.format(Title)
Profile_name = 'D:\\Master Astronomy Research year 2\\Master Project\\fig_profile_{}_Marta.png'.format(Title)

Fits_SNR_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Moment_zero_images\\{}_Marta_target.ms.split.contsub.cube.image.SNR_mom0.fits'.format(Title)
Fits_Ico_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Moment_zero_images\\{}_Marta_target.ms.split.contsub.cube.image.Ico_mom0.fits'.format(Title)
Fits_cube_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Profiles\\{}_Marta_target.ms.split.contsub.cube.image.fits'.format(Title)
Txt_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Profiles\\spectral_profile_{}_Marta.txt'.format(Title)

z = 4.341
D_L = 40340.3                # Mpc, use https://www.astro.ucla.edu/%7Ewright/CosmoCalc.html with H0 = 67.8, Om = 0.31 and Ovac = 0.69

Npts_5_SNR = 81
Npts_5_Ico = 81
Npts_12_Ico = 441
SNR_5 = 3.5566e-2
Ico_5 = 4.0083e-2
Ico_12 = 4.8961e-2
Abeam_SNR = 63.9442
Abeam_Ico = 63.8978


MakePlot(Fits_SNR_name, Npts_5_SNR, Abeam_SNR, Title, Cont_name, Savecondition)
MakeProfile(Txt_name, Fits_cube_name, Npts_5_SNR, Abeam_SNR, Title, Profile_name, Savecondition)
calc_param(Fits_SNR_name, Fits_Ico_name, Npts_5_SNR, Npts_5_Ico, Npts_12_Ico, SNR_5, Ico_5, Ico_12, Abeam_SNR, Abeam_Ico, z, D_L)


#%%

Title= 'AEG2'
Cont_name = 'D:\\Master Astronomy Research year 2\\Master Project\\fig_cont_{}_Marta.png'.format(Title)
Profile_name = 'D:\\Master Astronomy Research year 2\\Master Project\\fig_profile_{}_Marta.png'.format(Title)

Fits_SNR_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Moment_zero_images\\{}_Marta_target.ms.split.cube.image.SNR_mom0.fits'.format(Title)
Fits_Ico_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Moment_zero_images\\{}_Marta_target.ms.split.cube.image.Ico_mom0.fits'.format(Title)
Fits_cube_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Profiles\\{}_Marta_target.ms.split.cube.image.fits'.format(Title)
Txt_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Profiles\\spectral_profile_{}_Marta.txt'.format(Title)

z = 3.668
D_L = 33043.8                # Mpc, use https://www.astro.ucla.edu/%7Ewright/CosmoCalc.html with H0 = 67.8, Om = 0.31 and Ovac = 0.69

Npts_5_SNR = 79
Npts_5_Ico = 79
Npts_12_Ico = 239
SNR_5 = 8.2531e-2
Ico_5 = 1.0305e-1
Ico_12 = 1.2971e-1
Abeam_SNR = 39.8308
Abeam_Ico = 40.8144 


MakePlot(Fits_SNR_name, Npts_5_SNR, Abeam_SNR, Title, Cont_name, Savecondition)
MakeProfile(Txt_name, Fits_cube_name, Npts_5_SNR, Abeam_SNR, Title, Profile_name, Savecondition)
calc_param(Fits_SNR_name, Fits_Ico_name, Npts_5_SNR, Npts_5_Ico, Npts_12_Ico, SNR_5, Ico_5, Ico_12, Abeam_SNR, Abeam_Ico, z, D_L)




#%%



Title= 'AS2UDS026.0'
Cont_name = 'D:\\Master Astronomy Research year 2\\Master Project\\fig_cont_{}_Marta.png'.format(Title)
Profile_name = 'D:\\Master Astronomy Research year 2\\Master Project\\fig_profile_{}_Marta.png'.format(Title)

Fits_SNR_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Moment_zero_images\\{}_Marta_target.ms.split.cube.image.SNR_mom0.fits'.format(Title)
Fits_Ico_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Moment_zero_images\\{}_Marta_target.ms.split.cube.image.Ico_mom0.fits'.format(Title)
Fits_cube_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Profiles\\{}_Marta_target.ms.split.cube.image.fits'.format(Title)
Txt_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Profiles\\spectral_profile_{}_Marta.txt'.format(Title)

z = 3.296
D_L = 29088.7               # Mpc, use https://www.astro.ucla.edu/%7Ewright/CosmoCalc.html with H0 = 67.8, Om = 0.31 and Ovac = 0.69

Npts_5_SNR = 79
Npts_5_Ico = 79
Npts_12_Ico = 439
SNR_5 = 5.3054e-2
Ico_5 = 1.0333e-1
Ico_12 = 1.9051e-1
Abeam_SNR = 44.8977
Abeam_Ico = 44.9845


MakePlot(Fits_SNR_name, Npts_5_SNR, Abeam_SNR, Title, Cont_name, Savecondition)
MakeProfile(Txt_name, Fits_cube_name, Npts_5_SNR, Abeam_SNR, Title, Profile_name, Savecondition)
calc_param(Fits_SNR_name, Fits_Ico_name, Npts_5_SNR, Npts_5_Ico, Npts_12_Ico, SNR_5, Ico_5, Ico_12, Abeam_SNR, Abeam_Ico, z, D_L)


#%%


Title= 'CDFN1'
Cont_name = 'D:\\Master Astronomy Research year 2\\Master Project\\fig_cont_{}_NRAO.png'.format(Title)
Profile_name = 'D:\\Master Astronomy Research year 2\\Master Project\\fig_profile_{}_NRAO.png'.format(Title)

Fits_SNR_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Moment_zero_images\\{}_NRAO_target.ms.split.contsub.cube.image.SNR_mom0.fits'.format(Title)
Fits_Ico_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Moment_zero_images\\{}_NRAO_target.ms.split.contsub.cube.image.Ico_mom0.fits'.format(Title)
Fits_cube_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Profiles\\{}_NRAO_target.ms.split.contsub.cube.image.fits'.format(Title)
Txt_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Profiles\\spectral_profile_{}_NRAO.txt'.format(Title)

z = 3.149
D_L = 27543.4                       # Mpc, use https://www.astro.ucla.edu/%7Ewright/CosmoCalc.html with H0 = 67.8, Om = 0.31 and Ovac = 0.69

Npts_5_SNR = 79
Npts_5_Ico = 79
Npts_12_Ico = 439
SNR_5 = 1.5718e-1
Ico_5 = 1.5374e-1
Ico_12 = 1.6671e-1
Abeam_SNR = 31.1736
Abeam_Ico = 31.2147


MakePlot(Fits_SNR_name, Npts_5_SNR, Abeam_SNR, Title, Cont_name, Savecondition)
MakeProfile(Txt_name, Fits_cube_name, Npts_5_SNR, Abeam_SNR, Title, Profile_name, Savecondition)
calc_param(Fits_SNR_name, Fits_Ico_name, Npts_5_SNR, Npts_5_Ico, Npts_12_Ico, SNR_5, Ico_5, Ico_12, Abeam_SNR, Abeam_Ico, z, D_L)


#%%



Title= 'CDFN2'
Cont_name = 'D:\\Master Astronomy Research year 2\\Master Project\\fig_cont_{}_NRAO.png'.format(Title)
Profile_name = 'D:\\Master Astronomy Research year 2\\Master Project\\fig_profile_{}_NRAO.png'.format(Title)

Fits_SNR_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Moment_zero_images\\{}_NRAO_target.ms.split.contsub.cube.image.SNR_mom0.fits'.format(Title)
Fits_Ico_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Moment_zero_images\\{}_NRAO_target.ms.split.contsub.cube.image.Ico_mom0.fits'.format(Title)
Fits_cube_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Profiles\\{}_NRAO_target.ms.split.contsub.cube.image.fits'.format(Title)
Txt_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Profiles\\spectral_profile_{}_NRAO.txt'.format(Title)

z = 4.422
D_L = 41229.3                    # Mpc, use https://www.astro.ucla.edu/%7Ewright/CosmoCalc.html with H0 = 67.8, Om = 0.31 and Ovac = 0.69

Npts_5_SNR = 79
Npts_5_Ico = 79
Npts_12_Ico = 439
SNR_5 = 8.6521e-2
Ico_5 = 1.0196e-1
Ico_12 = 1.7195e-1
Abeam_SNR = 66.1418
Abeam_Ico = 66.2260


MakePlot(Fits_SNR_name, Npts_5_SNR, Abeam_SNR, Title, Cont_name, Savecondition)
MakeProfile(Txt_name, Fits_cube_name, Npts_5_SNR, Abeam_SNR, Title, Profile_name, Savecondition)
calc_param(Fits_SNR_name, Fits_Ico_name, Npts_5_SNR, Npts_5_Ico, Npts_12_Ico, SNR_5, Ico_5, Ico_12, Abeam_SNR, Abeam_Ico, z, D_L)




#%%

Title= 'AS2COS0031.1'
Cont_name = 'D:\\Master Astronomy Research year 2\\Master Project\\fig_cont_{}_NRAO.png'.format(Title)
Profile_name = 'D:\\Master Astronomy Research year 2\\Master Project\\fig_profile_{}_NRAO.png'.format(Title)

Fits_SNR_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Moment_zero_images\\{}_NRAO_target.ms.split.cube.image.SNR_mom0.fits'.format(Title)
Fits_Ico_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Moment_zero_images\\{}_NRAO_target.ms.split.cube.image.Ico_mom0.fits'.format(Title)
Fits_cube_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Profiles\\{}_NRAO_target.ms.split.cube.image.fits'.format(Title)
Txt_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Profiles\\spectral_profile_{}_NRAO.txt'.format(Title)

z = 3.643
D_L = 32776.1      # Mpc, use https://www.astro.ucla.edu/%7Ewright/CosmoCalc.html with H0 = 67.8, Om = 0.31 and Ovac = 0.69

Npts_5_SNR = 81
Npts_5_Ico = 81
Npts_12_Ico = 441
SNR_5 = 1.0527e-1
Ico_5 = 1.1789e-1
Ico_12 = 1.6867e-1
Abeam_SNR = 45.7982
Abeam_Ico = 45.8563


MakePlot(Fits_SNR_name, Npts_5_SNR, Abeam_SNR, Title, Cont_name, Savecondition)
MakeProfile(Txt_name, Fits_cube_name, Npts_5_SNR, Abeam_SNR, Title, Profile_name, Savecondition)
calc_param(Fits_SNR_name, Fits_Ico_name, Npts_5_SNR, Npts_5_Ico, Npts_12_Ico, SNR_5, Ico_5, Ico_12, Abeam_SNR, Abeam_Ico, z, D_L)






#%%

Title= 'AS2COS0054.1'
Cont_name = 'D:\\Master Astronomy Research year 2\\Master Project\\fig_cont_{}_NRAO.png'.format(Title)
Profile_name = 'D:\\Master Astronomy Research year 2\\Master Project\\fig_profile_{}_NRAO.png'.format(Title)

Fits_SNR_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Moment_zero_images\\{}_NRAO_target.ms.split.cube.image.SNR_mom0.fits'.format(Title)
Fits_Ico_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Moment_zero_images\\{}_NRAO_target.ms.split.cube.image.Ico_mom0.fits'.format(Title)
Fits_cube_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Profiles\\{}_NRAO_target.ms.split.cube.image.fits'.format(Title)
Txt_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Profiles\\spectral_profile_{}_NRAO.txt'.format(Title)

z = 3.174
D_L = 27805.5                         # Mpc, use https://www.astro.ucla.edu/%7Ewright/CosmoCalc.html with H0 = 67.8, Om = 0.31 and Ovac = 0.69

Npts_5_SNR = 81
Npts_5_Ico = 81
Npts_12_Ico = 441
SNR_5 = 1.1108e-1
Ico_5 = 1.3966e-1
Ico_12 = 1.3643e-1
Abeam_SNR = 46.4982
Abeam_Ico = 47.3124


MakePlot(Fits_SNR_name, Npts_5_SNR, Abeam_SNR, Title, Cont_name, Savecondition)
MakeProfile(Txt_name, Fits_cube_name, Npts_5_SNR, Abeam_SNR, Title, Profile_name, Savecondition)
calc_param(Fits_SNR_name, Fits_Ico_name, Npts_5_SNR, Npts_5_Ico, Npts_12_Ico, SNR_5, Ico_5, Ico_12, Abeam_SNR, Abeam_Ico, z, D_L)






























#%%


Title= 'AS2COS0001.1'
Cont_name = 'D:\\Master Astronomy Research year 2\\Master Project\\fig_cont_{}_NRAO.png'.format(Title)
Profile_name = 'D:\\Master Astronomy Research year 2\\Master Project\\fig_profile_{}_NRAO.png'.format(Title)

Fits_SNR_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Moment_zero_images\\{}_NRAO_combined_target.ms.split.cube.image.SNR_mom0.fits'.format(Title)
Fits_Ico_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Moment_zero_images\\{}_NRAO_combined_target.ms.split.cube.image.Ico_mom0.fits'.format(Title)
Fits_cube_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Profiles\\{}_NRAO_combined_target.ms.split.cube.image.fits'.format(Title)
Txt_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Profiles\\spectral_profile_{}_NRAO.txt'.format(Title)

z = 4.625
D_L = 43466.5      # Mpc, use https://www.astro.ucla.edu/%7Ewright/CosmoCalc.html with H0 = 67.8, Om = 0.31 and Ovac = 0.69

Npts_5_SNR = 81
Npts_5_Ico = 81
Npts_12_Ico = 441
SNR_5 = 6.3711e-2
Ico_5 = 9.3024e-2
Ico_12 = 1.6792e-1
Abeam_SNR = 69.6242
Abeam_Ico = 69.8101


MakePlot(Fits_SNR_name, Npts_5_SNR, Abeam_SNR, Title, Cont_name, Savecondition)
MakeProfile(Txt_name, Fits_cube_name, Npts_5_SNR, Abeam_SNR, Title, Profile_name, Savecondition)
calc_param(Fits_SNR_name, Fits_Ico_name, Npts_5_SNR, Npts_5_Ico, Npts_12_Ico, SNR_5, Ico_5, Ico_12, Abeam_SNR, Abeam_Ico, z, D_L)






#%%


Title= 'AS2COS0001.1_whole'
Cont_name = 'D:\\Master Astronomy Research year 2\\Master Project\\fig_cont_{}_NRAO.png'.format(Title)
Profile_name = 'D:\\Master Astronomy Research year 2\\Master Project\\fig_profile_{}_NRAO.png'.format(Title)

Fits_SNR_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Moment_zero_images\\as2cos1_sb_2_whole.cube.image.SNR_mom0.fits'.format(Title)
Fits_Ico_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Moment_zero_images\\as2cos1_sb_2_whole.cube.image.Ico_mom0.fits'.format(Title)
Fits_cube_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Profiles\\as2cos1_sb_2_whole.cube.image.fits'.format(Title)
Txt_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Profiles\\spectral_profile_as2cos1_sb_2_whole.txt'.format(Title)

z = 4.625
D_L = 43466.5      # Mpc, use https://www.astro.ucla.edu/%7Ewright/CosmoCalc.html with H0 = 67.8, Om = 0.31 and Ovac = 0.69

Npts_5_SNR = 81
Npts_5_Ico = 81
Npts_12_Ico = 441
SNR_5 = 6.4065e-2
Ico_5 = 1.0133e-1
Ico_12 = 2.2011e-1
Abeam_SNR = 69.4215
Abeam_Ico = 69.8148


MakePlot(Fits_SNR_name, Npts_5_SNR, Abeam_SNR, Title, Cont_name, Savecondition)
MakeProfile(Txt_name, Fits_cube_name, Npts_5_SNR, Abeam_SNR, Title, Profile_name, Savecondition)
calc_param(Fits_SNR_name, Fits_Ico_name, Npts_5_SNR, Npts_5_Ico, Npts_12_Ico, SNR_5, Ico_5, Ico_12, Abeam_SNR, Abeam_Ico, z, D_L)

#%%


Title= 'AS2COS0001.1_half'
Cont_name = 'D:\\Master Astronomy Research year 2\\Master Project\\fig_cont_{}_NRAO.png'.format(Title)
Profile_name = 'D:\\Master Astronomy Research year 2\\Master Project\\fig_profile_{}_NRAO.png'.format(Title)

Fits_SNR_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Moment_zero_images\\as2cos1_sb_2_half.cube.image.SNR_mom0.fits'.format(Title)
Fits_Ico_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Moment_zero_images\\as2cos1_sb_2_half.cube.image.Ico_mom0.fits'.format(Title)
Fits_cube_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Profiles\\as2cos1_sb_2_half.cube.image.fits'.format(Title)
Txt_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Profiles\\spectral_profile_as2cos1_sb_2_half.txt'.format(Title)

z = 4.625
D_L = 43466.5      # Mpc, use https://www.astro.ucla.edu/%7Ewright/CosmoCalc.html with H0 = 67.8, Om = 0.31 and Ovac = 0.69

Npts_5_SNR = 81
Npts_5_Ico = 81
Npts_12_Ico = 441
SNR_5 = 9.9907e-2
Ico_5 = 1.3081e-1
Ico_12 = 1.8728e-1
Abeam_SNR = 68.885
Abeam_Ico = 69.1465


MakePlot(Fits_SNR_name, Npts_5_SNR, Abeam_SNR, Title, Cont_name, Savecondition)
MakeProfile(Txt_name, Fits_cube_name, Npts_5_SNR, Abeam_SNR, Title, Profile_name, Savecondition)
calc_param(Fits_SNR_name, Fits_Ico_name, Npts_5_SNR, Npts_5_Ico, Npts_12_Ico, SNR_5, Ico_5, Ico_12, Abeam_SNR, Abeam_Ico, z, D_L)






















#%%

Fits_mom0_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Moment_zero_images\\as2uds10_128_05_100.split.cube.image.mom0.fits'
Fits_SNR_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Moment_zero_images\\as2uds10_128_05_100.split.cube.image.mom0.fits'
Fits_cube_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Profiles\\as2uds10_128_05_100.split.cube.image.fits'
Txt_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Profiles\\spectral_profile_as2uds10.txt'
Cont_name = 'D:\\Master Astronomy Research year 2\\Master Project\\fig_cont_AS2UDS010.0_v2.png'
Profile_name = 'D:\\Master Astronomy Research year 2\\Master Project\\fig_profile_AS2UDS010.0_v2.png'
Title= 'AS2UDS010.0'

z = 3.169
D_L = 27753.0       # Mpc, use https://www.astro.ucla.edu/%7Ewright/CosmoCalc.html with H0 = 67.8, Om = 0.31 and Ovac = 0.69

Npts_5_SNR = 79
Npts_5_Ico = 79
Npts_12_Ico = 439
SNR_5 = 0.13        #Jy km/s for 1 FWHM with radius 2.5
Ico_5 = 0.13        #Jy km/s for 2 FWM with radius 2.5
Ico_12 = 0.32       #Jy km/s for 2 FWHM with radius 6
Abeam_SNR = 54.7933
Abeam_Ico = 54.7933


MakePlot(Fits_SNR_name, Npts_5_SNR, Abeam_SNR, Title, Cont_name, Savecondition)
MakeProfile(Txt_name, Fits_cube_name, Npts_5_SNR, Abeam_SNR, Title, Profile_name, Savecondition)
calc_param(Fits_SNR_name, Fits_mom0_name, Npts_5_SNR, Npts_5_Ico, Npts_12_Ico, SNR_5, Ico_5, Ico_12, Abeam_SNR, Abeam_Ico, z, D_L)

#%%

Fits_mom0_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Moment_zero_images\\as2cos23-calib_128_05_100.split.cube.image.mom0.fits'
Fits_cube_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Profiles\\as2cos23-calib_128_05_100.split.cube.image.fits'
Txt_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Profiles\\spectral_profile_as2cos23.txt'
Cont_name = 'D:\\Master Astronomy Research year 2\\Master Project\\fig_cont_AS2COS0023.1.png'
Profile_name = 'D:\\Master Astronomy Research year 2\\Master Project\\fig_profile_AS2COS0023.1.png'
Title= 'AS2COS0023.1'

Npts = 1
Abeam = 1

RMS_AS2COS23 = MakePlot(Fits_mom0_name, Title, Cont_name, Savecondition)
FWHM_AS2COS23 = MakeProfile(Txt_name, Fits_cube_name, Title, Profile_name, Npts, Abeam, Savecondition)


#%%

Fits_mom0_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Moment_zero_images\\AS2COS54-my-calib.split.cube.image.mom0.fits'
Fits_cube_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Profiles\\AS2COS54-my-calib.split.cube.image.fits'
Txt_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Profiles\\spectral_profile_AS2COS54-my-calib.txt'
Cont_name = 'D:\\Master Astronomy Research year 2\\Master Project\\fig_cont_AS2COS54-my-calib.png'
Profile_name = 'D:\\Master Astronomy Research year 2\\Master Project\\fig_profile_AS2COS54-my-calib.png'
Title= 'AS2COS0054.1 my calib'

Npts = 81
Abeam = 53.4958

RMS_AS2COS54_my_calib = MakePlot(Fits_mom0_name, Title, Cont_name, Savecondition)
FWHM_AS2COS54_my_calib = MakeProfile(Txt_name, Fits_cube_name, Title, Profile_name, Npts, Abeam, Savecondition)

#%%

Fits_mom0_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Moment_zero_images\\AS2COS54-NRAO-calib.split.cube.image.Ico_mom0.fits'
Fits_SNR_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Moment_zero_images\AS2COS54-NRAO-calib.split.cube.image.SNR_mom0.fits'
Fits_cube_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Profiles\\AS2COS54-NRAO-calib.split.cube.image.fits'
Txt_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Profiles\\spectral_profile_AS2COS54-NRAO-calib.txt'
Cont_name = 'D:\\Master Astronomy Research year 2\\Master Project\\fig_cont_AS2COS54-NRAO-calib.png'
Profile_name = 'D:\\Master Astronomy Research year 2\\Master Project\\fig_profile_AS2COS54-NRAO-calib.png'
Title= 'AS2COS0054.1 NRAO calib'

z = 3.174
D_L = 27805.5       # Mpc, use https://www.astro.ucla.edu/%7Ewright/CosmoCalc.html with H0 = 67.8, Om = 0.31 and Ovac = 0.69

Npts_5_SNR = 81
Npts_5_Ico = 81
Npts_12_Ico = 441
SNR_5 = 8.425e-2       #Jy km/s for 1 FWHM with radius 2.5
Ico_5 = 1.1134e-1       #Jy km/s for 2 FWM with radius 2.5
Ico_12 = 1.0052e-1       #Jy km/s for 2 FWHM with radius 6
Abeam_SNR = 48.5645
Abeam_Ico = 49.4737


RMS_mom0 = MakePlot(Fits_mom0_name, Title, Cont_name, Savecondition)
FWHM = MakeProfile(Txt_name, Fits_cube_name, Title, Profile_name, Savecondition)
calc_param(Fits_SNR_name, Npts_5_SNR, Npts_5_Ico, Npts_12_Ico, SNR_5, Ico_5, Ico_12, Abeam_SNR, Abeam_Ico, RMS_mom0, z, D_L)

#%%

Fits_mom0_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Moment_zero_images\\AS2COS44-calib-only-NRAO.split.cube.image.mom0.fits'
Fits_cube_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Profiles\\AS2COS44-calib-only-NRAO.split.cube.image.fits'
Txt_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Profiles\\spectral_profile_AS2COS44_only_NRAO.txt'
Cont_name = 'D:\\Master Astronomy Research year 2\\Master Project\\fig_cont_AS2COS44_only_NRAO.png'
Profile_name = 'D:\\Master Astronomy Research year 2\\Master Project\\fig_profile_AS2COS44_only_NRAO.png'
Title= 'AS2COS0044.1 oN'

Npts = 81
Abeam = 44.891

RMS_AS2UD10 = MakePlot(Fits_mom0_name, Title, Cont_name, Savecondition)
FWHM_AS2UD10 = MakeProfile(Txt_name, Fits_cube_name, Title, Profile_name, Npts, Abeam, Savecondition)

#%%

Fits_mom0_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Moment_zero_images\\AS2UDS072.0_NRAO_target.split.cube.image.mom0.fits'
Fits_cube_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Profiles\\AS2UDS072.0_NRAO_target.split.cube.image.fits'
Txt_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Profiles\\spectral_profile_AS2UDS072.0_NRAO.txt'
Cont_name = 'D:\\Master Astronomy Research year 2\\Master Project\\fig_cont_AS2UDS072.0_NRAO.png'
Profile_name = 'D:\\Master Astronomy Research year 2\\Master Project\\fig_profile_AS2UDS072.0_NRAO.png'
Title= 'AS2UDS072.0 NRAO'

Npts = 79
Abeam = 33.7525

RMS_AS2UD10 = MakePlot(Fits_mom0_name, Title, Cont_name, Savecondition)
FWHM_AS2UD10 = MakeProfile(Txt_name, Fits_cube_name, Title, Profile_name, Npts, Abeam, Savecondition)


#%%


Fits_mom0_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Moment_zero_images\\AS2COS44-fully-my-calib.split.cube.image.mom0.fits'
Fits_cube_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Profiles\\AS2COS44-fully-my-calib.split.cube.image.fits'
Txt_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Profiles\\spectral_profile_AS2COS44_me.txt'
Cont_name = 'D:\\Master Astronomy Research year 2\\Master Project\\fig_cont_AS2COS44_me.png'
Profile_name = 'D:\\Master Astronomy Research year 2\\Master Project\\fig_profile_AS2COS44_me.png'
Title= 'AS2COS0044.1 me'

Npts = 79
Abeam = 46.5762

RMS_AS2UD10 = MakePlot(Fits_mom0_name, Title, Cont_name, Savecondition)
FWHM_AS2UD10 = MakeProfile(Txt_name, Fits_cube_name, Title, Profile_name, Npts, Abeam, Savecondition)
