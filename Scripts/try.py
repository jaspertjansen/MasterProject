# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 16:36:22 2023

@author: jansen
"""

#Check listobs

listobs(vis="as2uds10.ms")

# ------------------ Calibration --------------------------

# ------------------ Imaging --------------------------

# Step 1: Make a dirty image
# Expl: We use the dirty image to see if there is a signal and to get a 
#       feeling for the object we are about to study.

tclean(vis='as2uds10.ms',                       # Original file
       imagename='as2uds10_40_05.dirty',        # Destination file
       imsize=40,                               # How many pixels
       cell='0.5arcsec',                        # Size of pixel
       pblimit=-0.01,                           # Primary beam gain, - for no cor
       niter=0)                                 # No cleaning for dirty image

# Step 2: UV continuum subtraction
# Expl: We want to substract the continuum from the data. We first plot, so we
#       can find the line and hereafter substract the continuum, where we have
#       deselected the line. If no line is seen, select whole spectrum.

plotms(vis='as2uds10.ms',                       # Original file
       field='',                                # Field of source, after splitting 0
       ydatacolumn='data',                      # Data
       xaxis='channel',                         # Could also be freq
       yaxis='amp',                             # Amplitude
       correlation='RR',                        # ? Why only RR
       avgtime='1e8',                           # Some big timestep
       avgscan=True,                            # Averaging over all scans
       antenna='',                              # All antennas
       coloraxis='spw')                         # Spw's have color

# Optional: save spectrum in .png for freq and channel using the plotms viewer.

uvcontsub(vis='as2uds10.ms',                    # Original file
          fitspw='0~1:15~50',                   # Deselect line, skip if no line
          excludechans=True,                    # We want to exclude above chans
          want_cont=True)                       # Makes a continuum image

# Optional: save contsub and cont spectrum in .png for freq and channel using the plotms viewer.

# Step 3: Make a dirty cube of continuum-substracted data
# Expl: We want to see what has changed when we have substracted the continuum
#       and make a cube out of it.

tclean(vis='as2uds10.ms.contsub',               # Continuum substracted file
       imagename='as2uds10_40_05.contsub.cube_2MHz.',   # Destination file 
       imsize=40,                               # How many pixels 
       cell='0.5arcsec',                        # Size of pixel
       specmode='cube',                         # Cube for spectral line imaging
       pblimit=-0.01,                           # Primary beam gain, -for no cor 
       niter=0)                                 # No cleaning for dirty image
