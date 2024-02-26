# -*- coding: utf-8 -*-
"""
Created on Mon Feb  5 10:37:13 2024

@author: jansen
"""

# --------------------------------------------------------------------------- #
#                       Checking background sources
# --------------------------------------------------------------------------- #

# Step 1: Look for background source
# Expl: On dirty image look for background source in FoV of 3'x3'. If SNR>10,
#       we go on with the further process

tclean(vis='AS2COS0031.1_NRAO_target.ms.split',                 # Original file
       imagename='AS2COS0031.1_NRAO_target.ms.split.background2', # Destination file
       imsize=360,                              # 3'x3' FoV
       spw='2~13:5~58',                         # Remove noisyness
       cell='0.5arcsec',                        # Size of pixel
       pblimit=-0.01,                           # Primary beam gain, - for no cor
       niter=0)                                 # No cleaning for dirty image

# Step 2: Copy continuum MS file
# Expl: We do not want to mess up our data. Use in terminal:
#       cp -r AS2COS0031.1_NRAO_target.ms.split AS2COS0031.1_NRAO_target.ms.split.copy

# Step 3: Clean copied file
# Expl: Clean down to 1sigma? And save the model

tclean(vis='AS2COS0031.1_NRAO_target.ms.split.copy',               # Original file
       imagename='AS2COS0031.1_NRAO_target.ms.split.copy',     # Destination file
       imsize=360,                              # How many pixels
       cell='0.5arcsec',                        # Size of pixel
       specmode='mfs',                          # Line imaging
       weighting='uniform',                     # Best precision in location of source
       niter=1000000,                           # Some big number
       nsigma=1.0,                              # Stop cleaning at 1sigma
       fastnoise=False,                         # Noise estimation from data
       pblimit=-0.01,                           # Primary beam gain, - for no cor
       interactive=True,                        # We want to draw the mask ourself
       savemodel='modelcolumn')                 # We want to save the model into the MS file

# Step 4: uvsub
# Expl: subtract the model from the data to remove the background source

uvsub(vis='AS2COS0031.1_NRAO_target.ms.split.copy')

# Step 5: Split out corrected visibilities
# Expl: Split out the new visibilities in the corrected datacolumn

split(vis='AS2COS0031.1_NRAO_target.ms.split.copy',     # Original file
       outputvis='AS2COS0031.1_NRAO_target_nobckgrndsrc.ms.split',  # Destination file
       datacolumn='corrected')                  # Good data is in corrected column

# Step 6: Check background source
# Expl: Check if background source is gone

tclean(vis='AS2COS0031.1_NRAO_target_nobckgrndsrc.ms.split',                 # Original file
       imagename='AS2COS0031.1_NRAO_target_nobckgrndsrc.ms.split.background2', # Destination file
       imsize=360,                              # 3'x3' FoV
       spw='2~13:5~58',                         # Remove noisyness
       cell='0.5arcsec',                        # Size of pixel
       weighting='uniform',                     # Best precision in location of source
       pblimit=-0.01,                           # Primary beam gain, - for no cor
       niter=0)                                 # No cleaning for dirty image

# --------------------------------------------------------------------------- #
#                       uvplane fitting
# --------------------------------------------------------------------------- #

# Step 1: Single SPW
# Expl: We only want to have the SPW in which we expect the line

split(vis='AS2COS0028.1_NRAO_target.ms.split',  # Original file + _nobckgrndsrc for COS14, COS31 and CDFN8
      outputvis='AS2COS0028.1_NRAO_target.ms.split.line',   # Destinatio file
      spw=10,                                   # SPW where line is
      datacolumn='data')                        # Select data

# Step 2: Phaseshift
# Expl: In order to extremely certain that our low SNR data is messed up by
#       phase-errors, we shift the phasecenter to the center of emission. Check
#       for center in new dirty image.

tclean(vis='AS2COS0028.1_NRAO_target.ms.split.line',                 # Original file
       imagename='AS2COS0028.1_NRAO_target.ms.split.line.center', # Destination file
       imsize=128,                              # Smaller FoV
       cell='0.5arcsec',                        # Size of pixel
       pblimit=-0.01,                           # Primary beam gain, - for no cor
       niter=0)                                 # No cleaning for dirty image

myphasecenter='J2000 02h18m03.566s -04d55m27.214s'  # new phase-tracking center

phaseshift(vis='AS2COS0028.1_NRAO_target.ms.split.line',    # Original file
       outputvis='AS2COS0028.1_NRAO_target.ms.split.line.shifted',  # Destination file
       phasecenter=myphasecenter)               # Defined phasecenter

# Step 3: Run visbin
# Expl: First, specify binsize, then vis, then run the binning procedure and at
#       last the plotting and fitting functions.

uvbinsize = 5                                   # In Klambda

vis = 'AS2COS0028.1_NRAO_target.ms.split.line.shifted'  # Original file

execfile("/data1/jjansen/visbinning_2polarisations_new.py") # Binning procedure

execfile("/data1/jjansen/errorfit_new.py")          # Plotting and fitting procedure


















