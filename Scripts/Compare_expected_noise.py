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

File = 'D:\\Master Astronomy Research year 2\\Master Project\\Noise\\AS2COS44-calib-only-NRAO.split.noise.image.fits'

Data, Header = ReadFITS(File)

RMS = np.std(Data[0:52,0:128])
print ("RMS =", RMS*10**6, "uJy / beam")
