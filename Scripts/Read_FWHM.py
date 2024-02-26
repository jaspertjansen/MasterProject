# -*- coding: utf-8 -*-
"""
Created on Mon Oct 16 11:31:12 2023

@author: jansen
"""

from astropy.io import fits
import numpy as np

def ReadFITS(FileName):
    FITSFile = fits.open(FileName, lazy_load_hdu=True)
    Data = FITSFile[1].data
    Header = FITSFile[1].header
    FITSFile.close()
    return Data, Header

#%%
    
source_name = 'AS2UDS627.0'


Data, Header = ReadFITS('D:\\Master Astronomy Research year 2\\Master Project\\Table_with_FWHM2.fits')

index = np.where(source_name == Data['ID'])[0][0]

print('FWHM is:', Data['FWHM'][index], '+-', Data['FWHM_err'][index])

print('J_up', Data['J_up'][index])

print('Ico', Data['I_co'][index])

print('err Ico', Data['I_co_err'][index])

print('Lprime', Data['Lco_J'][index]/1e10)

print('err Lprime', Data['Lco_J_err'][index]/1e10)


print('z is:', Data['z_CO'][index])



