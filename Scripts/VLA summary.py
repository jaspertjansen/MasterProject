# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 16:36:22 2023

@author: jansen
"""



# ------------------ Calibration --------------------------

# ------------------ Imaging --------------------------

# Step 1: Check listobs
# Expl: Inspect the data and see if there is something weird or exceptional.

listobs(vis="as2uds10.ms")

# Step 2: Average over time
# Expl: Using the split function, we can average the data over time. Not much
#       data is lost and it greatly reduces computing times in the next steps!

split(vis='as2uds10.ms',                        # Original file
       outputvis='as2uds10.ms.split',           # Destination file
       datacolumn='data',                       # Original ms file may have data not in corrected column
       timebin='20s')                           # Time for which will be binned

# Step 3: Check listobs again
# Expl: Inspect the data and see if there is something weird or exceptional, again!
#       We can also see in which SPW the line is located using the redshift, as
#       restfreq=115.271203GHz/(1+z).

listobs(vis="as2uds10.ms.split")

# Step 4: Make a dirty continuum image
# Expl: We use the dirty image to see if there is a continuum by selecting the
#       spw's on both sides of the line spw. The image can be seen using imview.

tclean(vis='as2uds10.ms.split',                 # Original file
       imagename='as2uds10_128_05.split.dirty', # Destination file
       imsize=128,                              # How many pixels
       spw='9:5~58,11:5~58',                    # Two spw's, also mark out noice
       cell='0.5arcsec',                        # Size of pixel
       pblimit=-0.01,                           # Primary beam gain, - for no cor
       niter=0)                                 # No cleaning for dirty image

# Step 5: OPTIONAL: UV continuum subtraction
# Expl: OPTIONAL: We could substract the continuum from the data. Normally we would deselect
#       the line, but our lines are very faint. We deselected the outer spws and
#       the edge channels of spws, because these have much noise. Only do this is
#       the previous dirty image is not just noise.

uvcontsub(vis='as2uds10.ms.split',              # Original file
          fitspw='2~13:5~58',                   # Deselect noise at spw edges
          combine='spw',                        # As we have left out entire spw's, set combine
          excludechans=True,                    # We want to exclude above chans
          want_cont=False)                      # Makes a continuum image

# Step 6: Cleaning
# Expl: Using the tclean method in CASA, we can try to suppress the noise around
#       the oject. We draw a mask using the interactive cleaning method en clean
#       to a leven of 1sigma, because we have a low SNR.
#       IMPORTANT: SET GREEN CLEANING BAR TO CLEAN ALL CHANNELS AND POLARIZATIONS!

tclean(vis = 'as2uds10.ms.split',               # Original file (+.constsub if step 5 is taken)
       imagename='as2uds10_128_05_100.split.cube',   # Destination file
       imsize=128,                              # How many pixels
       cell='0.5arcsec',                        # Size of pixel
       width='100km/s',                         # Width of channel
       specmode='cube',                         # Line imaging
       spw='9~11',                              # spw's where line is located
       restfreq='27.649605GHz',                 # Set CO(1-0) restfreq=115.271203GHz/(1+z)
       reffreq='27.649605GHz',                  # Where 0km/s is defined, is the same
       mask='/data1/jjansen/MP_mask_6arcseconds.crtf',  # Predefined mask of 6" radius on centre pixel
       niter=1000000,                           # Some big number
       nsigma=1.0,                              # Stop cleaning at 1sigma
       fastnoise=False,                         # Noise estimation from data
       pblimit=-0.01,                           # Primary beam gain, - for no cor
       interactive=True)                        # We want to draw the mask ourself

# Step 7: Get out the profile
# Expl: Check for signal in the .cube.image file in imview spectral profile tool 
#       by drawing a region. Set x to velocity and y to flux density. Then, save
#       profile as .txt file. No function needed for this. Also, check which
#       channels have flux density for line. Needed in next step.

# Step 8: 0-th moment map
# Expl: We make an integrated intensity map by collapsing the cube. This should
#       look something like Marta's right figures in her figure 2 of Frias Castillo
#       et al. 2023.

immoments(imagename='as2uds10_128_05_100.split.cube.image',     # Original file
          chans='18~23',                        # Only collapse channels with signal
          moments=0,                            # Moment 0
          outfile='as2uds10_128_05_100.split.cube.image.mom0')  # Destination file

# Step 9: fits file
# Expl: Make a file to read to make a figure.

exportfits(imagename='as2uds10_128_05_100.split.cube.image.mom0',   # Original file
           fitsimage='as2uds10_128_05_100.split.cube.image.mom0.fits')  # Destination file














