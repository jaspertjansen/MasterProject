# -*- coding: utf-8 -*-
"""
Created on Tue Oct 10 15:30:16 2023

@author: jansen
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from matplotlib.ticker import LogFormatter
import matplotlib.patheffects as PathEffects
from matplotlib.patches import Ellipse  # this is for the beam contour
#plt.style.use('classic')

#%%

def ReadFITS(FileName):
    FITSFile = fits.open(FileName, lazy_load_hdu=True)
    Data = FITSFile[0].data[0][0]
    Header = FITSFile[0].header
    FITSFile.close()
    return Data, Header

def MakePlot(File, Title, Name):

    Offset = np.asarray((0,0,0,0))
    #read the fits files
    Data, Header = ReadFITS(File)
    Bmaj = Header['BMAJ']*3600.0    # major axis, convert deg -> arcsec
    Bmin = Header['BMIN']*3600.0
    Bpa = Header['BPA']+90

    print ("Beam:", Bmaj, Bmin, Bpa)

    PxScale = abs(Header['CDELT1']*3600.0) #arcsec

    DimX=128
    Extent0=np.asarray([-DimX*PxScale/2.,DimX*PxScale/2.,-DimX*PxScale/2.,DimX*PxScale/2.])+Offset
    #Extent0=np.asarray([-10,10,-10,10])+Offset

    print(-DimX*PxScale/2.)
    
    fig1 = plt.figure()
    ax=fig1.add_subplot(111)
    plt.xticks([])
    plt.yticks([])

    RMS = np.std(Data[0:52,0:128])
    print ("rms =", RMS)
    plt.xlim(-10,10)
    plt.ylim(-10,10)

    ax1=plt.imshow(Data*1000., cmap = 'viridis', origin='lower', interpolation = 'none', extent =Extent0)
    plt.contour(Data/RMS, levels=[-4,-2,2,4,6,8,10,12], colors =['0'], linewidths =[1], extent =Extent0)

    cb = plt.colorbar(ax1, fraction=0.035)
    cb.set_label('mJy/beam', fontsize = 12)
    # scale bar
    plt.text(8.25,-8.25, '2\"', color = '1',verticalalignment='center', horizontalalignment='center',backgroundcolor='none', fontsize = 16)
    plt.plot([7.25,9.25],[-9.25,-9.25], color = '1', lw = 2)
    #print(Bmaj, DimX*PxScale)
    Beam = Ellipse([-6.5, -6.5], width = Bmaj, height = Bmin, angle = Bpa, hatch = '/////', fc = 'none', ec = '1', lw=1, zorder = 10)
    ax.add_patch(Beam)

    #plt.text(-3.75, 3.65, Title, color = '1', verticalalignment='center', fontsize = 26)
    plt.text(-9.25, 8.75, Title, color = '1', verticalalignment='center', fontsize = 20)
    fig1.savefig(Name,bbox_inches='tight')
    plt.close(fig1)


#%%

File_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Moment_zero_images\\as2uds10_128_05_100.split.cube.image.mom0.fits'
Name = 'D:\\Master Astronomy Research year 2\\Master Project\\fig_cont_AS2UDS010.0.png'
Title= 'AS2UDS010.0'

MakePlot(File_name, Title, Name)

#%%

File_name = 'D:\\Master Astronomy Research year 2\\Master Project\\Moment_zero_images\\as2cos23-calib_128_05_100.split.cube.image.mom0.fits'
Name = 'D:\\Master Astronomy Research year 2\\Master Project\\fig_cont_AS2COS0023.1.png'
Title= 'AS2COS0023.1'

MakePlot(File_name, Title, Name)











