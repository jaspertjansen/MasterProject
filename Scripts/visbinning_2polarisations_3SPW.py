# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 10:36:10 2024

@author: jansen
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 10:23:41 2024

@author: jansen
"""

# visbinning_cyc1.py
# Written by James Simpson (originally called visbinning_jms.py)
# Downloaded by JH on 01 Dec 2015
# Modified by Matus 2018

# Matus: get the UV-plane bin size (in klambda) from the CASA environment
# specified via a variable uvbinsize (all small, to avoid typecase errors)
# check if variable initialized, and larger than 0:
if (uvbinsize >= 0.0):
    print ("UV bin size = ", uvbinsize, "klambda")
else:
    print ("Warning: need to specify uvbinsize")
    exit()


import numpy as np				# JH added, as was missing
import matplotlib.pyplot as plt			# JH added, as was missing
from matplotlib.backends.backend_pdf import PdfPages
import os.path
import subprocess
from scipy.optimize import curve_fit		# JH added

################# Define code for grabbing visibilities ###################
###########################################################################
def get_vis(inms=None):

    #speed of light
    c  = 2.99792458E8

    #Open measurement set
    

    #Initialize empty arrays

    vreal =  np.zeros(0)
    vimag =  np.zeros(0)
    uvdist = np.zeros(0)
    weight = np.zeros(0)


    #Loop over each SPW
    No_SPWs = 3
    for j in range(No_SPWs):
        ms.open(inms)

        print ("MS file loaded")

        # Load data for SPW j
        ms.selectinit(datadescid=j)
        d = ms.getdata(['u','v','real','imaginary','uvdist','axis_info','flag','weight'],ifraxis=True)
        
        ms.close()


        # Sloppy coding, but make array for uvdistances matched in size to visibility data
        uvdist_temp = np.zeros(np.shape(d['real']))
        weight_temp = np.zeros(np.shape(d['real']))

        # Get the frequency of each channel
        freq = d['axis_info']['freq_axis']['chan_freq'][:]
        conv_fac = (c/freq)       # Conversion to wavelength units

        # Now fill temporary array (0 - XX pol, 1 - YY pol) and convert to klambda
        print ("Here", np.shape(freq)[0])
        for i in range(np.shape(freq)[0]):                # Loop through frequency array
            uvdist_temp[0,i,:,:] = d['uvdist']/conv_fac[i]   # Convert uvdistance to lambda (XX pol)
            # !!!!! ONLY 1 polarisation in NOEMA !!!!
            #uvdist_temp[0,i,:,:] = d['uvdist']/conv_fac[i]   # Convert uvdistance to lambda (YY pol)
            uvdist_temp[1,i,:,:] = d['uvdist']/conv_fac[i]   # Convert uvdistance to lambda (YY pol)
            weight_temp[:,i,:,:] = d['weight']

        # Note: vreal, vimag are given in Stokes XX, YY, _NOT_ Stokes I!
        mask = np.where(d['flag'] != True)              # Identify data that has NOT been flagged
        vreal = np.append( vreal,  (d['real'][mask]).flatten() )        # Extract good data, and collapse multi-D arrays  - Real
        vimag = np.append( vimag,  (d['imaginary'][mask]).flatten() )   # Extract good data, and collapse multi-D arrays  - Imaginary
        uvdist = np.append( uvdist, ((uvdist_temp[mask]).flatten()) )   # Extract good data, and collapse multi-D arrays  - uv distance
        weight = np.append( weight, ((weight_temp[mask]).flatten()) )

    # use more "human" units
    uvdist = uvdist*1.E-3      # Convert lambda to klambda
    vreal = vreal*1.E3      # convert Jy to mJy
    vimag = vimag*1.E3

    # print some statistics. Useful for debugging
    #print 'Mean of real: '+ str(np.mean(vreal))+'mJy'
    #print 'Mean of imaginary: '+ str(np.mean(vimag))+'mJy'
    #print 'Median uv distance: '+ "{0:.2f}".format(np.median(uvdist))+'klambda (Min: '+"{0:.2f}".format(min(uvdist))+', Max: '+"{0:.2f}".format(max(uvdist))+')'



    return(uvdist,vreal,vimag,weight)

################# Define code for binning visibilities ####################
###########################################################################
def bin_vis(vis=None, uvdis=None, binmin= None, binmax= None, binwidth= None, weight=None):
#Returns the mean visibility and mean uvdistance with respect to the user specified bins

    #Create bin limits (set from data if not specified)
    if binmin is None:  binmin = min(uvdis)
    if binmax is None:  binmax = max(uvdis)

    nbins = int((binmax-binmin)/binwidth)+1     # Determine the number of bins required

    # Setup empty arrays
    vis_bin = np.zeros(nbins,dtype=float)
    vis_meanerr_bin = np.zeros(nbins,dtype=float)
    vis_sigma_bin = np.zeros(nbins,dtype=float)
    uv_bin  = np.zeros(nbins,dtype=float)
    vis_num = np.zeros(nbins,dtype=int)

    # Loop through data and bin visibilities
    for i in range(0, nbins):
        x = np.where(np.logical_and(uvdist > binmin+(i*binwidth), uvdist < binmin+((i+1)*binwidth) ) )      # Indices of vis at correct uv dist


        vis_bin[i] = np.sum(vis[x]*weight[x])/(np.sum(weight[x]))         # Mean of vis at correct uv dist
        tmp = vis[x]
        tmp = np.sort(tmp)
        vis_meanerr_bin[i] = np.std(vis[x])/(np.shape(x)[1])**0.5   # Standard error on mean
        vis_sigma_bin[i] = np.std(vis[x])    # Store 1-sigma scatter in visibilities
        uv_bin[i] = np.mean(uvdis[x])        # Mean uvdist
        vis_num[i] = np.shape(x)[1]          # Number of visibilities in the bin
    return(uv_bin, vis_bin, vis_meanerr_bin, vis_sigma_bin, vis_num )



####################################################################################
######################### Start of main code #######################################
####################################################################################

######################### Define initial parameters #################################


#bin width in klambda
# grabbed from the CASA interface
binwidth_uv = uvbinsize


################################## END ############################################


##################### Read in split MS files ####################
#Open split_MS_phaseshifted directory and read in file names
msfiles=[]
msfiles = np.asarray(vis)
########################### END ###################################


################ Start loop to bin visibilities ###################


# Subroutine that returns the uv distance, visibilities and weights
#uvdist,vreal,vimag,weights = get_vis(inms=msfiles[i].split('/')[-1])
uvdist,vreal,vimag,weights = get_vis(inms=vis)

# Call subroutines to bin the real and imaginary visibilities
a,b,c,d,e = bin_vis(vis=vreal,uvdis=uvdist,weight=weights,binmin=min(uvdist),binmax=max(uvdist),binwidth=binwidth_uv)
f,g,h,k,l = bin_vis(vis=vimag,uvdis=uvdist,weight=weights,binmin=min(uvdist),binmax=max(uvdist),binwidth=binwidth_uv)

print ("Writing into a file")
# Write the binned amp+errors to a text file for future analysis
output='#UVdist (kl) '+'\t'+'Real Amp '+'\t'+'Error on Real mean '+'\t'+'Imag Amp '+'\t'+'Error on Imag mean '+'\t'+'No of binned visibilities \n'

for j in range(0,np.shape(a)[0]):
    output+= str(a[j])+'\t'+str(b[j])+'\t'+str(c[j])+'\t'+str(g[j])+'\t'+str(h[j])+'\t'+str(e[j])+'\n'   

outfile = open(vis+'_binned_'+(str)(uvbinsize)+'_klambda.txt' ,'w')
outfile.write(output)
outfile.close()

########################### END of vis analysis #######################
