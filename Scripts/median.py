# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 12:08:21 2024

@author: jansen
"""

import numpy as np

def calc_median(x,x_err):
    
    list_of_medians = []   
    
    for j in range(10000):
        new_x = []
        for i in range(len(x)):
            new_x.append(np.random.normal(x[i],x_err[i]))
        list_of_medians.append(np.nanmedian(new_x))
        
    return np.nanmedian(list_of_medians), np.nanstd(list_of_medians)


#%%

name = ["AS2COS0001.1", "AS2COS0002.1", "AS2COS0006.1", "AS2COS0008.1", "AS2COS0009.1", "AS2COS0011.1", "AS2COS0013.1", "AS2COS0014.1", \
        "AS2COS0023.1", "AS2COS0028.1", "AS2COS0031.1", "AS2COS0044.1", "AS2COS0054.1", "AS2COS0065.1", "AS2COS0066.1", "AS2COS0139.1", \
        "AS2UDS009.0", "AS2UDS010.0", "AS2UDS011.0", "AS2UDS012.0", "AS2UDS014.0", "AS2UDS026.0", "AS2UDS072.0", "AS2UDS126.0", \
        "AS2UDS231.0", "AEG2", "AEG3", "CDFN1", "CDFN2", "CDFN8"]

Ico_1_0 = [0.24, np.NaN, np.NaN, 0.19, np.NaN, np.NaN, 0.44, 0.23, 0.07, 0.15, 0.17, np.NaN, \
       0.14, 0.10, 0.18, 0.18, 0.12, 0.34, 0.11, 0.22, np.NaN, 0.19, np.NaN, 0.22, 0.12, 0.13, np.NaN, 0.17, 0.17, 0.17]
Ico_1_0_err = [0.15, np.NaN, np.NaN, 0.10, np.NaN, np.NaN, 0.08, 0.09, 0.05, 0.10, 0.06, np.NaN, \
           0.11, 0.15, 0.09, 0.09, 0.08, 0.09, 0.05, 0.07, np.NaN, 0.07, np.NaN, 0.12, 0.09, 0.06, np.NaN, 0.14, 0.07, 0.07]
       
Lco_1_0 = [19.8, np.NaN, np.NaN, 10.4, np.NaN, np.NaN, 14.3, 9.1, 5.4, 6.6, 9.5, np.NaN, 6.2, 3.0, 8.5, \
           8.8, 4.9, 15.5, 7.1, 6.7, np.NaN, 9.0, np.NaN, 6.4, 5.2, 7.4, np.NaN, 7.5, 13.2, 11.8]
Lco_1_0_err = [12.3, np.NaN, np.NaN, 5.3, np.NaN, np.NaN, 2.8, 3.4, 3.8, 4.5, 3.1, np.NaN, 5.0, 4.4, 4.2, \
             4.3, 3.2, 4.0, 3.5, 2.1, np.NaN, 3.4, np.NaN, 3.3, 4.0, 3.5, np.NaN, 6.3, 5.3, 4.0]

r_31 = [np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, 0.91, 0.68, np.NaN, 0.68, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, \
        1.08, 0.67, np.NaN, 0.59, np.NaN, np.NaN, np.NaN, 0.67, np.NaN, np.NaN, np.NaN, 0.72, np.NaN, np.NaN]
r_31_err = [np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, 0.18, 0.27, np.NaN, 0.48, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, \
            0.75, 0.19, np.NaN, 0.22, np.NaN, np.NaN, np.NaN, 0.37, np.NaN, np.NaN, np.NaN, 0.62, np.NaN, np.NaN]
        
r_41 = [np.NaN, np.NaN, np.NaN, 0.69, np.NaN, np.NaN, np.NaN, np.NaN, 1.69, np.NaN, 0.68, np.NaN, 0.71, 2.74, 0.56, 0.89, np.NaN, \
        np.NaN, 0.55, np.NaN, np.NaN, 0.48, np.NaN, np.NaN, 0.75, 0.34, np.NaN, np.NaN, 0.22, 0.36]
r_41_err = [np.NaN, np.NaN, np.NaN, 0.36, np.NaN, np.NaN, np.NaN, np.NaN, 1.20, np.NaN, 0.23, np.NaN, 0.60, 4.04, 0.30, 0.44, np.NaN, \
            np.NaN, 0.30, np.NaN, np.NaN, 0.2, np.NaN, np.NaN, 0.70, 0.16, np.NaN, np.NaN, 0.09, 0.26]
        
r_51 = [0.14, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, \
        np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN]
r_51_err = [0.09, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, \
            np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN]

cont = [np.NaN, 30.29, np.NaN, np.NaN, np.NaN, 22.96, np.NaN, np.NaN, 17.70, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, \
        np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, 14.42, np.NaN, np.NaN, np.NaN, np.NaN, 31.58, 29.89, np.NaN]

z = [4.625, 4.595, 4.62, 3.581, 2.26, 4.786, 2.608, 2.921, 4.341, 3.097, 3.643, 2.58, 3.174, 2.414, 3.247, 3.292, 2.942, \
     3.169, 4.073, 2.52, 3.804, 3.296, 2.406, 2.436, 3.119, 3.668, 4.032, 3.149, 4.422, 4.144]

FWHM = [1311, np.NaN, 347, 771, np.NaN, np.NaN, 210, 574, 85, 605, 430, np.NaN, 726, 294, 222, 439, 477, 483, 836, \
        862, np.NaN, 203, np.NaN, 462, 395, 2415, np.NaN, 502, 475, 857]
FWHM_err = [563, np.NaN, 231, 183, np.NaN, np.NaN, 21, 100, 1626, 145, 79, np.NaN, 221, 188, 73, 98, 227, 97, 339, \
            197, np.NaN, 167, np.NaN, 115, 229, 671, np.NaN, 135, 80, 338]

SNR = [1.9, 1.2, 0.8, 4.1, 0.2, 0.4, 8.5, 6.1, 1.9, 4.6, 6.5, 0.2, 3.2, 1.5, 2.3, 3.9, 1.8, 4.3, 1.9, 4.6, -0.2, 1.7, \
       1.1, 4.2, 1.8, 3.1, -0.1, 3.9, 4.3, 2.3]

l = [3, np.NaN, np.NaN, 11, np.NaN, np.NaN, 8, 6, 26, 6, 11, np.NaN, 11, 24, 9, 14, 10, 6, 9, 5, np.NaN, 8, np.NaN, 6, 9, 5, np.NaN, 7, 3, 6]
l_err = [2, np.NaN, np.NaN, 6, np.NaN, np.NaN, 2, 2, 19, 4, 4, np.NaN, 9, 36, 5, 7, 7, 2, 5, 2, np.NaN, 3, np.NaN, 4, 9, 2, np.NaN, 6, 1, 4]

RMS_cube = [143, 141, 132, 85, 92, 101, 81, 56, 65, 77, 50, 116, 82, 170, 107, 88, 120, 97, 71, \
            54, 81, 88, 89, 86, 122, 53, 41, 87, 57, 65]

noise = [41.77, 46.70, 30.78, 25.94, np.NaN, 27.89, np.NaN, 14.67, np.NaN, 27.95, 18.70, 33.87, 23.46, 51.67, 58.06, 26.29, \
         28.03, 28.29, 19.68, 16.89, np.NaN, 24.04, 38.32, 25.06, 39.59, 13.30, np.NaN, 30.97, 22.28, 19.47]

#%%


print("Cont = ", np.nanmedian(cont), "+-", np.nanstd(cont))

print("z = ", np.nanmedian(z), "+-", np.nanstd(z), "low", np.nanmin(z), "high", np.nanmax(z))

FWHM_mean, FWHM_std = calc_median(FWHM, FWHM_err)
print("FWHM = ", FWHM_mean, "+-", FWHM_std)

print("SNR = ", np.nanmedian(SNR), "+-", np.nanstd(SNR))
 
Ico_1_0_mean, Ico_1_0_std = calc_median(Ico_1_0, Ico_1_0_err)
print("Ico_1_0 = ", Ico_1_0_mean, "+-", Ico_1_0_std)

l_mean, l_std = calc_median(l, l_err)
print("l = ", l_mean, "+-", l_std)

Lco_1_0_mean, Lco_1_0_std = calc_median(Lco_1_0, Lco_1_0_err)
print("Lco_1_0 = ", Lco_1_0_mean, "+-", Lco_1_0_std)

r_31_mean, r_31_std = calc_median(r_31, r_31_err)
print("r_31 = ", r_31_mean, "+-", r_31_std)

r_41_mean, r_41_std = calc_median(r_41, r_41_err)
print("r_41 = ", r_41_mean, "+-", r_41_std)

r_51_mean, r_51_std = calc_median(r_51, r_51_err)
print("r_51 = ", r_51_mean, "+-", r_51_std)

print("RMS cube = ", np.nanmedian(RMS_cube), "+-", np.nanstd(RMS_cube))

print("noise = ", np.nanmedian(noise), "+-", np.nanstd(noise))






#%%
f = open("Jasperdata.txt", "w")

for i in range(len(name)):
    f.write("{} {} {}\n".format(name[i], z[i], Ico[i]))
#    f.write()


#%%
    
np.genfromtxt("Jasperdata.txt", delimiter="", dtype=None)


#%%
labels = np.genfromtxt('Jasperdata.txt', delimiter='', usecols=0, dtype=str)
raw_data = np.genfromtxt('Jasperdata.txt', delimiter='')[:,1:]

print(labels)
print(raw_data[np.where(labels=="AS2COS0001.1")[0][0]])
print(raw_data[:,1])


np.where(labels=="AS2COS0001.1")[0][0]






#%%

name = ["AS2COS0001.1", "AS2COS0002.1", "AS2COS0006.1", "AS2COS0008.1", "AS2COS0009.1", "AS2COS0011.1", "AS2COS0013.1", "AS2COS0014.1", \
        "AS2COS0023.1", "AS2COS0028.1", "AS2COS0031.1", "AS2COS0044.1", "AS2COS0054.1", "AS2COS0065.1", "AS2COS0066.1", "AS2COS0139.1", \
        "AS2UDS009.0", "AS2UDS010.0", "AS2UDS011.0", "AS2UDS012.0", "AS2UDS014.0", "AS2UDS026.0", "AS2UDS072.0", "AS2UDS126.0", \
        "AS2UDS231.0", "AEG2", "AEG3", "CDFN1", "CDFN2", "CDFN8"]

Ico_1_0 = [0.24, 0.06, 0.23, 0.19, np.NaN, 0.08, 0.44, 0.23, 0.07, 0.15, 0.17, 0.08, \
       0.14, 0.10, 0.18, 0.18, 0.12, 0.34, 0.11, 0.22, np.NaN, 0.19, 0.03, 0.22, 0.12, 0.13, np.NaN, 0.17, 0.17, 0.17]
Ico_1_0_err = [0.15, 0.09, 0.11, 0.10, np.NaN, 0.13, 0.08, 0.09, 0.05, 0.10, 0.06, 0.14, \
           0.11, 0.15, 0.09, 0.09, 0.08, 0.09, 0.05, 0.07, np.NaN, 0.07, 0.06, 0.12, 0.09, 0.06, np.NaN, 0.14, 0.07, 0.07]
       
Lco_1_0 = [19.8, 4.6, 18.6, 10.4, np.NaN, 7.1, 14.3, 9.1, 5.4, 6.6, 9.5, 2.6, 6.2, 3.0, 8.5, \
           8.8, 4.9, 15.5, 7.1, 6.7, np.NaN, 9.0, 1.0, 6.4, 5.2, 7.4, np.NaN, 7.5, 13.2, 11.8]
Lco_1_0_err = [12.3, 6.9, 8.9, 5.3, np.NaN, 11.7, 2.8, 3.4, 3.8, 4.5, 3.1, 4.4, 5.0, 4.4, 4.2, \
             4.3, 3.2, 4.0, 3.5, 2.1, np.NaN, 3.4, 1.8, 3.3, 4.0, 3.5, np.NaN, 6.3, 5.3, 4.0]

r_31 = [np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, 0.91, 0.68, np.NaN, 0.68, np.NaN, 1.94, np.NaN, np.NaN, np.NaN, np.NaN, \
        1.08, 0.67, np.NaN, 0.59, np.NaN, np.NaN, 4.56, 0.67, np.NaN, np.NaN, np.NaN, 0.72, np.NaN, np.NaN]
r_31_err = [np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, 0.18, 0.27, np.NaN, 0.48, np.NaN, 3.35, np.NaN, np.NaN, np.NaN, np.NaN, \
            0.75, 0.19, np.NaN, 0.22, np.NaN, np.NaN, 8.46, 0.37, np.NaN, np.NaN, np.NaN, 0.62, np.NaN, np.NaN]
        
r_41 = [np.NaN, np.NaN, np.NaN, 0.69, np.NaN, np.NaN, np.NaN, np.NaN, 1.69, np.NaN, 0.68, np.NaN, 0.71, 2.74, 0.56, 0.89, np.NaN, \
        np.NaN, 0.55, np.NaN, np.NaN, 0.48, np.NaN, np.NaN, 0.75, 0.34, np.NaN, np.NaN, 0.22, 0.36]
r_41_err = [np.NaN, np.NaN, np.NaN, 0.36, np.NaN, np.NaN, np.NaN, np.NaN, 1.20, np.NaN, 0.23, np.NaN, 0.60, 4.04, 0.30, 0.44, np.NaN, \
            np.NaN, 0.30, np.NaN, np.NaN, 0.2, np.NaN, np.NaN, 0.70, 0.16, np.NaN, np.NaN, 0.09, 0.26]
        
r_51 = [0.14, 0.17, 0.42, np.NaN, np.NaN, 0.81, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, \
        np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN]
r_51_err = [0.09, 0.31, 0.21, np.NaN, np.NaN, 1.33, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, \
            np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN]

cont = [np.NaN, 30.29, np.NaN, np.NaN, np.NaN, 22.96, np.NaN, np.NaN, 17.70, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, \
        np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, 14.42, np.NaN, np.NaN, np.NaN, np.NaN, 31.58, 29.89, np.NaN]

z = [4.625, 4.595, 4.62, 3.581, 2.26, 4.786, 2.608, 2.921, 4.341, 3.097, 3.643, 2.58, 3.174, 2.414, 3.247, 3.292, 2.942, \
     3.169, 4.073, 2.52, 3.804, 3.296, 2.406, 2.436, 3.119, 3.668, 4.032, 3.149, 4.422, 4.144]

FWHM = [1311, np.NaN, 347, 771, np.NaN, np.NaN, 210, 574, 85, 605, 430, np.NaN, 726, 294, 222, 439, 477, 483, 836, \
        862, np.NaN, 203, np.NaN, 462, 395, 2415, np.NaN, 502, 475, 857]
FWHM_err = [563, np.NaN, 231, 183, np.NaN, np.NaN, 21, 100, 1626, 145, 79, np.NaN, 221, 188, 73, 98, 227, 97, 339, \
            197, np.NaN, 167, np.NaN, 115, 229, 671, np.NaN, 135, 80, 338]

SNR = [1.9, 1.2, 0.8, 4.1, 0.2, 0.4, 8.5, 6.1, 1.9, 4.6, 6.5, 0.2, 3.2, 1.5, 2.3, 3.9, 1.8, 4.3, 1.9, 4.6, -0.2, 1.7, \
       1.1, 4.2, 1.8, 3.1, -0.1, 3.9, 4.3, 2.3]

l = [3, 21, 11, 11, np.NaN, 21, 8, 6, 26, 6, 11, 17, 11, 24, 9, 14, 10, 6, 9, 5, np.NaN, 8, 40, 6, 9, 5, np.NaN, 7, 3, 6]
l_err = [2, 32, 5, 6, np.NaN, 34, 2, 2, 19, 4, 4, 30, 9, 36, 5, 7, 7, 2, 5, 2, np.NaN, 3, 75, 4, 9, 2, np.NaN, 6, 1, 4]

RMS_cube = [143, 141, 132, 85, 92, 101, 81, 56, 65, 77, 50, 116, 82, 170, 107, 88, 120, 97, 71, \
            54, 81, 88, 89, 86, 122, 53, 41, 87, 57, 65]

noise = [41.77, 46.70, 30.78, 25.94, np.NaN, 27.89, np.NaN, 14.67, np.NaN, 27.95, 18.70, 33.87, 23.46, 51.67, 58.06, 26.29, \
         28.03, 28.29, 19.68, 16.89, np.NaN, 24.04, 38.32, 25.06, 39.59, 13.30, np.NaN, 30.97, 22.28, 19.47]

