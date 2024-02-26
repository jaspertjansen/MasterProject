# -*- coding: utf-8 -*-
"""
Created on Mon Nov  6 20:48:32 2023

@author: jansen
"""

from astropy.io import fits
import numpy as np


def ReadFITS(FileName):
    FITSFile = fits.open(FileName, lazy_load_hdu=True)
    Data = FITSFile[0].data[0][0]
    Header = FITSFile[0].header
    FITSFile.close()
    return Data, Header

File = 'D:\\Master Astronomy Research year 2\\Master Project\\Noise\\CDFN8_NRAO_target.ms.split.noise.fits'
#File = 'D:\\Master Astronomy Research year 2\\Master Project\\Noise\\as2cos1_sb_2_whole.noise.image.fits'
#File = 'D:\\Master Astronomy Research year 2\\Master Project\\Noise\\CDFN1_NRAO_target.ms.split.noise.fits'



freq_bw = 52.1003     # Found in proposal in MHz
exp_noise = 36.1596      # Found in proposal in uJy/beam

Data, Header = ReadFITS(File)

RMS_54_chan = np.std(Data[0:52,0:128])
#RMS_54_chan = np.std(Data[76:128,0:128])
RMS_freq_chan = np.sqrt(54/(freq_bw/2)) * RMS_54_chan
print("RMS =", RMS_freq_chan*10**6, "uJy / beam")
print("Fraction of exp noise", RMS_freq_chan*10**6/exp_noise)
