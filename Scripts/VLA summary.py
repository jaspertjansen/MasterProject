# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 16:36:22 2023

@author: jansen
"""



# ------------------ Calibration by myself --------------------------
# Done for AS2COS54

# Step 1: Download from Archive
# Expl: Search for '21A-254' or Hodge on the VLA archive (https://data.nrao.edu/portal/#/)
#       and select the dataset (using naming file from Marta). Set the calibration 
#       to uncalibrated!!! Also note the observing datae, is needed in step 3.

# Step 2: Unpack the .tar file.
# Expl: Unpack the tar file outside of CASA in a folder with the name of the source
#       using to get the .ms file: tar zxvf AS2COS54-uncalib.tar.gz and 
#       tar zxvf 21A-254.sb39393798.eb39561560.59308.094199803236.ms.tgz


# Step 3: Observing logs
# Expl: Check the observing logs using http://www.vla.nrao.edu/cgi-bin/oplogs.cgi
#       by selecting the observing date and inspect the log and keep this in mind
#       for later.

# Step 4: Check listobs
# Expl: Inspect the data using the listobs function in CASA, note down which
#       source is which calibrator

listobs(vis="21A-254.sb39393798.eb39561560.59308.094199803236.ms")

# Step 5: Antennas
# Expl: Plot the antennas. We use this to choose a reference antenna. Choose an
#       antenna close to centre, but not too close that it suffers from shadowing.
#       Choose antenna ea05, unless it appears in the log in step 3.

# Step 6: Inspect elevation
# Expl: We check the elevation over time of the target and flux calibrator. If
#       they are off, we have to correct in the future. We will be strongly dependant
#       on the opactity and gaincurve corrections.

plotms(vis='21A-254.sb39393798.eb39561560.59308.094199803236.ms', # Original file
       xaxis='time',                            # x-axis
       yaxis='elevation',                       # y-axis
       correlation='RR,LL',                     # All polarizations
       avgchannel='64',                         # 64 channels in a spw
       spw='0:4~60',                            # Just choose a random spw
       coloraxis='field')                       # Sources have different color

# Step 7: Inspect amplitude
# Expl: Zoom in and look for 0 values. If these values are found, make a box and
#       check them in the logger. There could be a problem with one of the antennas.

plotms(vis='21A-254.sb39393798.eb39561560.59308.094199803236.ms',   # Original file
       xaxis='time',                            # x-axis
       yaxis='amp',                             # y-axis
       correlation='RR,LL',                     # All polarizations
       avgchannel='64',                         # 64 channels in a spw
       spw='0:4~60',                            # Just choose a random spw
       coloraxis='field')                       # Sources have different color

# Step 8: OPTIONAL: Flag antennas
# Expl: OPTIONAL: If misbehaving antennas appear in the data, flag them.

flagdata(vis='21A-254.sb39393798.eb39561560.59308.094199803236.ms',     # Original file
         mode='list',                           # Can flag multiple
         inpfile=["field='2,3' antenna='ea12' timerange='03:41:00~04:10:00'",   # All spw's
                  "field='2,3' antenna='ea07,ea08' timerange='03:21:40~04:10:00' spw='1'"]) # Multiple antennas

# Step 9: Update antenna positions
# Expl: Update the antenna positions using online access.

gencal(vis='21A-254.sb39393798.eb39561560.59308.094199803236.ms',   # Original file
       caltable='antpos.cal',                   # Table to write to
       caltype='antpos',                        # Antenna positions
       antenna='')                              # All antennas

# Step 10: Antenna efficiencies
# Expl: The gaincurve describes how each antenna behaves as a function of elevation, 
#       for each receiver band.

gencal(vis='21A-254.sb39393798.eb39561560.59308.094199803236.ms',   # Original file
       caltable='gaincurve.cal',                # Table to write to
       caltype='gceff')                         # Efficiency

# Step 11: Opacities
# Expl: Only needed for high frequencies, we mostly have the Ka-band, so needed.
#       Atmospheric opacity affects flux density scale for different elevations.
#       First, generate weather conditions, then make a table.

myTau = plotweather(vis='21A-254.sb39393798.eb39561560.59308.094199803236.ms',  # Original file
                    doPlot=True)                # Return a plot image file

gencal(vis='21A-254.sb39393798.eb39561560.59308.094199803236.ms',   # Original file
       caltable='opacity.cal',                  # Table to write to
       caltype='opac',                          # Opacity
       spw='0~17',                              # All spw
       parameter=myTau)                         # Write to this parameter

# Step 12: Model for flux calibrator
# Expl: We will scale the observed flux calibrator with the model of the flux
#       calibrator. We use the Ka (A) model, as we observe in the Ka band. This
#       could depend on which source/redshift we are observing.

setjy(vis='21A-254.sb39393798.eb39561560.59308.094199803236.ms',   # Original file,
      listmodels=True)                          # Look at all the models

setjy(vis='21A-254.sb39393798.eb39561560.59308.094199803236.ms',   # Original file,
      field='2',                                # Field ID of flux calibrator
      spw='0~15',                               # For all spw's
      scalebychan=True,                         # For channels different scaling
      model='3C286_A.im')                       # 3C286 or 3C48 in Ka (A) band

# Step 13: Phase variations
# Expl: Inspect if we need to correct the phase variations with time before solving 
#       for the bandpass to prevent decorrelation of the bandpass solution. Most 
#       of the times, the phase variation as a function of channel is modest, so 
#       we can average over all channels to increase the SNR.

plotms(vis='21A-254.sb39393798.eb39561560.59308.094199803236.ms',   # Original file
       field='2',                               # Field ID of bandpass calibrator
       xaxis='time',                            # x-axis
       yaxis='phase',                           # y-axis
       correlation='RR',                        # Polarization
       avgchannel='64',                         # Average over channels
       antenna='ea05&ea26')                     # Refant and random

# Step 14: Delay calibration
# Expl: The delay is the linear slope of phase across frequency. They were almost
#       horizontal, but it is still a good idea to do it.

gaincal(vis='21A-254.sb39393798.eb39561560.59308.094199803236.ms',  # Original file
        caltable='delays.cal',                  # Table to write to
        field='2',                              # Field ID of bandpass calibrator
        refant='ea05',                          # Refant
        gaintype='K',                           # Delay
        gaintable=['antpos.cal','gaincurve.cal','opacity.cal'])     # Opacity might not be available

# Step 15: Bandpass calibration
# Expl: First correct for decorrelation, and inspect the result, then apply phase 
#       solution.

gaincal(vis='21A-254.sb39393798.eb39561560.59308.094199803236.ms',  # Original file
        caltable='bpphase.gcal',                # Table to write to
        field='2',                              # Field ID of bandpass calibrator
        spw='0~15:20~40',                       # All spw, no edges
        refant='ea05',                          # Refant
        calmode='p',                            # Correct for decorrelation
        solint='int',                           # Integration time of 10s
        minsnr=2.0,                             # Only good solutions
        gaintable=['antpos.cal','gaincurve.cal','delays.cal', 'opacity.cal'])  # Opacity might not be available

plotms(vis='bpphase.gcal',                      # Solution to inspect
       gridrows=3,                              # Figure rows
       gridcols=3,                              # Figure columns
       xaxis='time',                            # x-axis
       yaxis='phase',                           # y-axis
       iteraxis='antenna',                      # Subplot iteration
       coloraxis='corr',                        # RR, LL, etc.
       plotrange=[0,0,-180,180])                # Full phase

bandpass(vis='21A-254.sb39393798.eb39561560.59308.094199803236.ms',     # Original dfile
         caltable='bandpass.bcal',              # Table to write to
         field='2',                             # Field ID of bandpass calibrator
         refant='ea05',                         # Refant
         solint='inf',                          # Average over scan
         solnorm=False,                         # Do not normalize
         gaintable=['antpos.cal','gaincurve.cal','delays.cal','bpphase.gcal', 'opacity.cal'])   # Opacity might not be available

# Step 16: Phase solution on integrated time
# Expl: Solve and inspect for antenna-based phase and amplitude gain calibration.

gaincal(vis='21A-254.sb39393798.eb39561560.59308.094199803236.ms',  # Original file
        caltable='intphase.gcal',               # Table to write to
        field='0,2',                            # Amplitude/gain and bandpass calibrator
        spw='0~15:4~60',                        # All spws, no edges
        solint='int',                           # Integrated time
        refant='ea05',                          # Refant
        minsnr=2.0,                             # Only good solutions
        calmode='p',                            # Phase
        gaintable=['antpos.cal','gaincurve.cal','delays.cal','bandpass.bcal', 'opacity.cal'])   # Opacity might not be available

plotms(vis='intphase.gcal',                     # Solution to inspect
       gridrows=3,                              # Figure rows
       gridcols=3,                              # Figure columns
       xaxis='time',                            # x-axis
       yaxis='phase',                           # y-axis
       iteraxis='antenna',                      # Subplot iteration
       coloraxis='corr',                        # RR, LL, etc.
       plotrange=[0,0,-180,180])                # Full phase

# Step 17: Phase solution on a scan time
# Expl: To be used for the target source. And inspect solution.

gaincal(vis='21A-254.sb39393798.eb39561560.59308.094199803236.ms',  # Original file
        caltable='scanphase.gcal',              # Table to write to
        field='0,2',                            # Amplitude/gain and bandpass calibrator
        spw='0~15:4~60',                        # All spws, no edges
        solint='inf',                           # Scan time
        refant='ea05',                          # Refant
        minsnr=2.0,                             # Only good solutions
        calmode='p',                            # Phase
        gaintable=['antpos.cal','gaincurve.cal','delays.cal','bandpass.bcal', 'opacity.cal'])   # Opacity might not be available

plotms(vis='scanphase.gcal',                    # Solution to inspect
       gridrows=3,                              # Figure rows
       gridcols=3,                              # Figure columns
       xaxis='time',                            # x-axis
       yaxis='phase',                           # y-axis
       iteraxis='antenna',                      # Subplot iteration
       coloraxis='corr',                        # RR, LL, etc.
       plotrange=[0,0,-180,180])                # Full phase

# Step 18: Amplitude solutions
# Expl: Apply bandpass and step 16 and look at results. Plot will show the residual
#       phase error. If scatter is large, need to flag more data.

gaincal(vis='21A-254.sb39393798.eb39561560.59308.094199803236.ms',  # Original file
        caltable='amp.gcal',                    # Table to write to
        field='0,2',                            # Amplitude/gain and bandpass calibrator
        spw='0~15:4~60',                        # All spws, no edges
        solint='inf',                           # Scan time
        refant='ea05',                          # Refant
        minsnr=2.0,                             # Only good solutions
        calmode='ap',                           # Amplitude
        gaintable=['antpos.cal','gaincurve.cal','opacity.cal','delays.cal','bandpass.bcal','intphase.gcal'])    # Opacity might not be available

plotms(vis='amp.gcal',                          # Solution to inspect
       gridrows=3,                              # Figure rows
       gridcols=3,                              # Figure columns
       xaxis='time',                            # x-axis
       yaxis='phase',                           # y-axis
       iteraxis='antenna',                      # Subplot iteration
       coloraxis='corr',                        # RR, LL, etc.
       plotrange=[-1,-1,-20,20])                # Range of residual phase

# Step 19: Flux calibration
# Expl: Use flux calibrator to derive the flux of the other calibrators.

fluxscale(vis='21A-254.sb39393798.eb39561560.59308.094199803236.ms',    # Original file
          caltable='amp.gcal',                  # Table to write to
          fluxtable='flux.cal',                 # Table to write to
          reference='2',                        # Field ID of flux calibrator
          incremental=True)                     # With increments

# Step 20: Apply calibration
# Expl: First for gain/phase calibrator, then for flux and bandpass calibrator,
#       then the target source. Opacity has been omitted.

applycal(vis='21A-254.sb39393798.eb39561560.59308.094199803236.ms',     # Original file
         field='0',                             # Gain/phase calibrator
         gaintable=['antpos.cal','gaincurve.cal','opacity.cal','delays.cal','bandpass.bcal','intphase.gcal','amp.gcal','flux.cal'],
         gainfield=['','','','2','2','0','0','0'], # Field ID's for tables 
         calwt=False)                           # Do not weigh?

applycal(vis='21A-254.sb39393798.eb39561560.59308.094199803236.ms',     # Original file
         field='2',                             # Flux/bandpass calibrator
         gaintable=['antpos.cal','gaincurve.cal','opacity.cal','delays.cal','bandpass.bcal','intphase.gcal','amp.gcal','flux.cal'],   # Opacity omitted
         gainfield=['','','','2','2','2','2','2'], # Field ID's for tables
         calwt=False)                           # Do not weigh?

applycal(vis='21A-254.sb39393798.eb39561560.59308.094199803236.ms',     # Original file
         field='1',                             # Target source field
         gaintable=['antpos.cal','gaincurve.cal','opacity.cal','delays.cal','bandpass.bcal','scanphase.gcal','amp.gcal','flux.cal'],  # Opacity omitted
         gainfield=['','','','2','2','0','0','0'], # Field ID's for tables
         calwt=False)                           # Do not weigh?

# Step 21: Inspect calibration
# Expl: Look at corrected amp vs time and search for weird signals

plotms(vis='21A-254.sb39393798.eb39561560.59308.094199803236.ms',   # Original file
       xaxis='time',                            # x-axis
       yaxis='amp',                             # y-axis
       ydatacolumn='corrected',                 # Take new data column
       field='2',                               # Field ID of bandpass calibrator
       spw='10:4~60',                           # Look at spw 10
       correlation='RR,LL',                     # Polarizations
       avgchannel='64',                         # Average all channels in spw
       coloraxis='antenna1')                    # Color by antenna

plotms(vis='21A-254.sb39393798.eb39561560.59308.094199803236.ms',   # Original file
       xaxis='uvdist',                          # x-axis
       yaxis='amp',                             # y-axis
       ydatacolumn='corrected',                 # Take new data column
       field='1',                               # Field ID of target
       spw='10:4~60',                           # Look at spw 10
       correlation='RR,LL',                     # Polarizations
       avgchannel='64',                         # Average all channels in spw
       #antenna='!ea07;!ea12;!ea23', 
       coloraxis='antenna2')                    # Color by antenna

# Step 22: OPTIONAL: Flag data
# Expl: OPTIONAL: Can be done manually in plotms for every spw. Can also be done
#       directly in CASA with the following command

flagdata(vis='21A-254.sb39393798.eb39561560.59308.094199803236.ms',     # Original file 
         mode='list',                           # List is most easy
         inpfile=["antenna='ea25'"])            # Flag an antenna for the whole observation

# Step 23: OPTIONAL: Rerun tables
# Expl: OPTIONAL: After flagging, the calibration must be reapplied. 

gaincal(vis='21A-254.sb39393798.eb39561560.59308.094199803236.ms',  # Original file
        caltable='bpphase.gcal',                # Table to write to
        field='2',                              # Field ID of bandpass calibrator
        spw='0~15:20~40',                       # All spw, no edges
        refant='ea05',                          # Refant
        calmode='p',                            # Correct for decorrelation
        solint='int',                           # Integration time of 10s
        minsnr=2.0,                             # Only good solutions
        gaintable=['antpos.cal','gaincurve.cal','delays.cal', 'opacity.cal'])  # Opacity might not be available

bandpass(vis='21A-254.sb39393798.eb39561560.59308.094199803236.ms',     # Original dfile
         caltable='bandpass.bcal',              # Table to write to
         field='2',                             # Field ID of bandpass calibrator
         refant='ea05',                         # Refant
         solint='inf',                          # Average over scan
         solnorm=False,                         # Do not normalize
         gaintable=['antpos.cal','gaincurve.cal','delays.cal','bpphase.gcal', 'opacity.cal'])   # Opacity might not be available

gaincal(vis='21A-254.sb39393798.eb39561560.59308.094199803236.ms',  # Original file
        caltable='intphase.gcal',               # Table to write to
        field='0,2',                            # Amplitude/gain and bandpass calibrator
        spw='0~15:4~60',                        # All spws, no edges
        solint='int',                           # Integrated time
        refant='ea05',                          # Refant
        minsnr=2.0,                             # Only good solutions
        calmode='p',                            # Phase
        gaintable=['antpos.cal','gaincurve.cal','delays.cal','bandpass.bcal', 'opacity.cal'])   # Opacity might not be available

gaincal(vis='21A-254.sb39393798.eb39561560.59308.094199803236.ms',  # Original file
        caltable='scanphase.gcal',              # Table to write to
        field='0,2',                            # Amplitude/gain and bandpass calibrator
        spw='0~15:4~60',                        # All spws, no edges
        solint='inf',                           # Scan time
        refant='ea05',                          # Refant
        minsnr=2.0,                             # Only good solutions
        calmode='p',                            # Phase
        gaintable=['antpos.cal','gaincurve.cal','delays.cal','bandpass.bcal', 'opacity.cal'])   # Opacity might not be available

gaincal(vis='21A-254.sb39393798.eb39561560.59308.094199803236.ms',  # Original file
        caltable='amp.gcal',                    # Table to write to
        field='0,2',                            # Amplitude/gain and bandpass calibrator
        spw='0~15:4~60',                        # All spws, no edges
        solint='inf',                           # Scan time
        refant='ea05',                          # Refant
        minsnr=2.0,                             # Only good solutions
        calmode='ap',                           # Amplitude
        gaintable=['antpos.cal','gaincurve.cal','opacity.cal','delays.cal','bandpass.bcal','intphase.gcal'])    # Opacity might not be available

fluxscale(vis='21A-254.sb39393798.eb39561560.59308.094199803236.ms',    # Original file
          caltable='amp.gcal',                  # Table to write to
          fluxtable='flux.cal',                 # Table to write to
          reference='2',                        # Field ID of flux calibrator
          incremental=True)                     # With increments

applycal(vis='21A-254.sb39393798.eb39561560.59308.094199803236.ms',     # Original file
         field='0',                             # Gain/phase calibrator
         gaintable=['antpos.cal','gaincurve.cal','opacity.cal','delays.cal','bandpass.bcal','intphase.gcal','amp.gcal','flux.cal'],
         gainfield=['','','','2','2','0','0','0'], # Field ID's for tables 
         calwt=False)                           # Do not weigh?

applycal(vis='21A-254.sb39393798.eb39561560.59308.094199803236.ms',     # Original file
         field='2',                             # Flux/bandpass calibrator
         gaintable=['antpos.cal','gaincurve.cal','opacity.cal','delays.cal','bandpass.bcal','intphase.gcal','amp.gcal','flux.cal'],   # Opacity omitted
         gainfield=['','','','2','2','2','2','2'], # Field ID's for tables
         calwt=False)                           # Do not weigh?

applycal(vis='21A-254.sb39393798.eb39561560.59308.094199803236.ms',     # Original file
         field='1',                             # Target source field
         gaintable=['antpos.cal','gaincurve.cal','opacity.cal','delays.cal','bandpass.bcal','scanphase.gcal','amp.gcal','flux.cal'],  # Opacity omitted
         gainfield=['','','','2','2','0','0','0'], # Field ID's for tables
         calwt=False)                           # Do not weigh?

# Step 24: Split out the object
# Expl: Only take the data of the target object

split(vis='21A-254.sb39393798.eb39561560.59308.094199803236.ms',    # Original file
      outputvis='AS2COS54-my-calib.ms',         # New file name
      field='1',                                # Field ID of target
      spw='0~15')                               # Only science spws









# ------------------ Calibration by NRAO --------------------------
# Done for AS2COS54

# Step 1: Download from Archive
# Expl: Search for '21A-254' or Hodge on the VLA archive (https://data.nrao.edu/portal/#/)
#       and select the dataset (using naming file from Marta). Set the calibration 
#       to calibrated!!! 

# Step 2: Unpack the .tar file.
# Expl: Unpack the tar file outside of CASA in a folder with the name of the source
#       using to get the .ms file: tar zxvf 21A-254.sb39658500.eb39706045.59359.146427581014.tar.gz 
#       and tar zxvf 21A-254.sb39658500.eb39706045.59359.146427581014.ms.tgz

# Step 3: Check listobs
# Expl: Inspect the data using the listobs function in CASA, note down which
#       source is which calibrator, check observing dates and RA and Dec to
#       check if source is the correct source

listobs(vis="21A-254.sb39658500.eb39706045.59359.146427581014.ms")

# Step 4: Inspect calibration
# Expl: Look at corrected amp vs time and uvdist and search for weird signals and 
#       loop trough all fields and spw's

plotms(vis='21A-254.sb39658500.eb39706045.59359.146427581014.ms',   # Original file
       xaxis='time',                            # x-axis
       yaxis='amp',                             # y-axis
       ydatacolumn='corrected',                 # Take new data column
       field='0',                               # Field ID of bandpass calibrator
       spw='0',                                 # Look at spw 0
       correlation='RR,LL',                     # Polarizations
       avgchannel='64',                         # Average all channels in spw
       coloraxis='antenna1')                    # Color by antenna

plotms(vis='21A-254.sb39658500.eb39706045.59359.146427581014.ms',   # Original file
       xaxis='uvdist',                          # x-axis
       yaxis='amp',                             # y-axis
       ydatacolumn='corrected',                 # Take new data column
       field='0',                               # Field ID of target
       spw='0',                                 # Look at spw 0
       correlation='RR,LL',                     # Polarizations
       avgchannel='64',                         # Average all channels in spw
       #antenna='!ea07;!ea12;!ea23', 
       coloraxis='antenna2')                    # Color by antenna

# Step 5: OPTIONAL: Flag data
# Expl: OPTIONAL: Can be done manually in plotms for every spw. Can also be done
#       directly in CASA with the following command

flagdata(vis='21A-254.sb39658500.eb39706045.59359.146427581014.ms',
         mode='list', 
         inpfile=["field='0' antenna='ea01,ea21,ea18'",
                  "field='0' antenna='ea11' spw='0'",
                  "field='0' antenna='ea27' spw='8,9,10'"])

# Step 6: OPTIONAL: Rerun tables
# Expl: OPTIONAL: After flagging, the calibration must be reapplied. Antpos is
#       not available for NRAO calibrated data.

gaincal(vis='21A-254.sb39658500.eb39706045.59359.146427581014.ms',  # Original file
        caltable='bpphase.gcal',                # Table to write to
        field='2',                              # Field ID of bandpass calibrator
        spw='0~15:20~40',                       # All spw, no edges
        refant='ea05',                          # Refant
        calmode='p',                            # Correct for decorrelation
        solint='int',                           # Integration time of 10s
        minsnr=2.0,                             # Only good solutions
        gaintable=['unknown.session_1.caltables/21A-254.sb39658500.eb39706045.59359.146427581014.ms.hifv_priorcals.s5_2.gc.tbl',            # Gaincurve
                   'unknown.session_1.caltables/21A-254.sb39658500.eb39706045.59359.146427581014.ms.hifv_finalcals.s13_2.finaldelay.tbl',   # Delays
                   'unknown.session_1.caltables/21A-254.sb39658500.eb39706045.59359.146427581014.ms.hifv_priorcals.s5_3.opac.tbl'])         # Opacity

bandpass(vis='21A-254.sb39658500.eb39706045.59359.146427581014.ms',     # Original dfile
         caltable='bandpass.bcal',              # Table to write to
         field='2',                             # Field ID of bandpass calibrator
         refant='ea05',                         # Refant
         solint='inf',                          # Average over scan
         solnorm=False,                         # Do not normalize
         gaintable=['unknown.session_1.caltables/21A-254.sb39658500.eb39706045.59359.146427581014.ms.hifv_priorcals.s5_2.gc.tbl',            # Gaincurve
                    'unknown.session_1.caltables/21A-254.sb39658500.eb39706045.59359.146427581014.ms.hifv_finalcals.s13_2.finaldelay.tbl',   # Delays
                    'bpphase.gcal', 
                    'unknown.session_1.caltables/21A-254.sb39658500.eb39706045.59359.146427581014.ms.hifv_priorcals.s5_3.opac.tbl'])         # Opacity

gaincal(vis='21A-254.sb39658500.eb39706045.59359.146427581014.ms',  # Original file
        caltable='intphase.gcal',               # Table to write to
        field='0,2',                            # Amplitude/gain and bandpass calibrator
        spw='0~15',                             # All spws, no edges
        solint='int',                           # Integrated time
        refant='ea05',                          # Refant
        minsnr=2.0,                             # Only good solutions
        calmode='p',                            # Phase
        gaintable=['unknown.session_1.caltables/21A-254.sb39658500.eb39706045.59359.146427581014.ms.hifv_priorcals.s5_2.gc.tbl',            # Gaincurve
                   'unknown.session_1.caltables/21A-254.sb39658500.eb39706045.59359.146427581014.ms.hifv_finalcals.s13_2.finaldelay.tbl',   # Delays
                   'bpphase.gcal', 
                   'unknown.session_1.caltables/21A-254.sb39658500.eb39706045.59359.146427581014.ms.hifv_priorcals.s5_3.opac.tbl'])         # Opacity

gaincal(vis='21A-254.sb39658500.eb39706045.59359.146427581014.ms',  # Original file
        caltable='scanphase.gcal',              # Table to write to
        field='0,2',                            # Amplitude/gain and bandpass calibrator
        spw='0~15',                             # All spws
        solint='inf',                           # Scan time
        refant='ea05',                          # Refant
        minsnr=2.0,                             # Only good solutions
        calmode='p',                            # Phase
        gaintable=['unknown.session_1.caltables/21A-254.sb39658500.eb39706045.59359.146427581014.ms.hifv_priorcals.s5_2.gc.tbl',            # Gaincurve
                   'unknown.session_1.caltables/21A-254.sb39658500.eb39706045.59359.146427581014.ms.hifv_finalcals.s13_2.finaldelay.tbl',   # Delays
                   'bpphase.gcal', 
                   'unknown.session_1.caltables/21A-254.sb39658500.eb39706045.59359.146427581014.ms.hifv_priorcals.s5_3.opac.tbl'])         # Opacity

gaincal(vis='21A-254.sb39658500.eb39706045.59359.146427581014.ms',  # Original file
        caltable='amp.gcal',                    # Table to write to
        field='0,2',                            # Amplitude/gain and bandpass calibrator
        spw='0~15',                             # All spws
        solint='inf',                           # Scan time
        refant='ea05',                          # Refant
        minsnr=2.0,                             # Only good solutions
        calmode='ap',                           # Amplitude
        gaintable=['unknown.session_1.caltables/21A-254.sb39658500.eb39706045.59359.146427581014.ms.hifv_priorcals.s5_2.gc.tbl',            # Gaincurve
                   'unknown.session_1.caltables/21A-254.sb39658500.eb39706045.59359.146427581014.ms.hifv_priorcals.s5_3.opac.tbl',         # Opacity
                   'unknown.session_1.caltables/21A-254.sb39658500.eb39706045.59359.146427581014.ms.hifv_finalcals.s13_2.finaldelay.tbl',   # Delays
                   'bpphase.gcal',
                   'intphase.gcal'])   

fluxscale(vis='21A-254.sb39658500.eb39706045.59359.146427581014.ms',    # Original file
          caltable='amp.gcal',                  # Table to write to
          fluxtable='flux.cal',                 # Table to write to
          reference='2',                        # Field ID of flux calibrator
          incremental=True)                     # With increments

applycal(vis='21A-254.sb39658500.eb39706045.59359.146427581014.ms',     # Original file
         field='0',                             # Gain/phase calibrator
         gaintable=['unknown.session_1.caltables/21A-254.sb39658500.eb39706045.59359.146427581014.ms.hifv_priorcals.s5_2.gc.tbl',            # Gaincurve
                    'unknown.session_1.caltables/21A-254.sb39658500.eb39706045.59359.146427581014.ms.hifv_priorcals.s5_3.opac.tbl',         # Opacity
                    'unknown.session_1.caltables/21A-254.sb39658500.eb39706045.59359.146427581014.ms.hifv_finalcals.s13_2.finaldelay.tbl',   # Delays
                    'bandpass.bcal',
                    'intphase.gcal',
                    'amp.gcal',
                    'flux.cal'],
         gainfield=['','','2','2','0','0','0'], # Field ID's for tables 
         calwt=False)                           # Do not weigh?

applycal(vis='21A-254.sb39658500.eb39706045.59359.146427581014.ms',     # Original file
         field='2',                             # Flux/bandpass calibrator
         gaintable=['unknown.session_1.caltables/21A-254.sb39658500.eb39706045.59359.146427581014.ms.hifv_priorcals.s5_2.gc.tbl',            # Gaincurve
                    'unknown.session_1.caltables/21A-254.sb39658500.eb39706045.59359.146427581014.ms.hifv_priorcals.s5_3.opac.tbl',          # Opacity
                    'unknown.session_1.caltables/21A-254.sb39658500.eb39706045.59359.146427581014.ms.hifv_finalcals.s13_2.finaldelay.tbl',   # Delays
                    'bandpass.bcal',
                    'intphase.gcal',
                    'amp.gcal',                    
                    'flux.cal'],
         gainfield=['','','2','2','2','2','2'], # Field ID's for tables
         calwt=False)                           # Do not weigh?

applycal(vis='21A-254.sb39658500.eb39706045.59359.146427581014.ms',     # Original file
         field='1',                             # Target source field
         gaintable=['unknown.session_1.caltables/21A-254.sb39658500.eb39706045.59359.146427581014.ms.hifv_priorcals.s5_2.gc.tbl',            # Gaincurve
                    'unknown.session_1.caltables/21A-254.sb39658500.eb39706045.59359.146427581014.ms.hifv_priorcals.s5_3.opac.tbl',          # Opacity
                    'unknown.session_1.caltables/21A-254.sb39658500.eb39706045.59359.146427581014.ms.hifv_finalcals.s13_2.finaldelay.tbl',   # Delays
                    'bandpass.bcal',
                    'scanphase.gcal',
                    'amp.gcal',
                    'flux.cal'],  # Opacity omitted
         gainfield=['','','2','2','0','0','0'], # Field ID's for tables
         calwt=False)                           # Do not weigh?

# Step 7: Split out the object
# Expl: Only take the data of the target object

split(vis='21A-254.sb39658500.eb39706045.59359.146427581014.ms',    # Original file
      outputvis='AS2COS44-calib-NRAO-me.ms',    # New file name
      field='1',                                # Field ID of target
      spw='0~15')                               # Only science spws









# ------------------ Imaging --------------------------
# Done for AS2UDS10.ms, received from Marta. Both observations were already
# calibrated and merged.

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

# Step 4: Compare expected noise with the observed noise
# Expl: The expected noise is stated in the observing proposal. If there is a big
#       difference, the NRAO calibration is off, or something else is going wrong

tclean(vis = 'as2uds10.ms.split',               # Original file (+.constsub if step 5 is taken)
       imagename='as2uds10_128_05_100.split.noise',     # Destination file
       imsize=128,                              # How many pixels
       cell='0.5arcsec',                        # Size of pixel
       width='100km/s',                         # Width of channel
       specmode='mfs',                         # Line imaging
       spw='9:5~58,11:5~58',                              # spw's where line is located
       restfreq='87.5386MHz',                 # Set CO(1-0) restfreq=115.271203GHz/(1+z)
       reffreq='87.5386MHz',                  # Where 0km/s is defined, is the same
       mask='/data1/jjansen/MP_mask_6arcseconds.crtf',  # Predefined mask of 6" radius on centre pixel
       niter=1000000,                           # Some big number
       nsigma=1.0,                              # Stop cleaning at 1sigma
       fastnoise=False,                         # Noise estimation from data
       pblimit=-0.01,                           # Primary beam gain, - for no cor
       interactive=True)                        # We want to draw the mask ourself


# Step 5: Export and compare
# Expl: Make fits file and load this into "Compare_expected_noise.py"

exportfits(imagename='AS2COS44-calib-NRAO-me.split.noise.image',   # Original file 
           fitsimage='AS2COS44-calib-NRAO-me.split.noise.image.fits')  # Destination file 

# Step 6: Make a dirty continuum image
# Expl: We use the dirty image to see if there is a continuum by selecting the
#       spw's on both sides of the line spw. The image can be seen using imview.

tclean(vis='as2uds10.ms.split',                 # Original file
       imagename='as2uds10_128_05.split.dirty', # Destination file
       imsize=128,                              # How many pixels
       spw='9:5~58,11:5~58',                    # Two spw's, also mark out noice
       cell='0.5arcsec',                        # Size of pixel
       pblimit=-0.01,                           # Primary beam gain, - for no cor
       niter=0)                                 # No cleaning for dirty image

# Step 7: OPTIONAL: UV continuum subtraction
# Expl: OPTIONAL: We could substract the continuum from the data. Normally we would deselect
#       the line, but our lines are very faint. We deselected the outer spws and
#       the edge channels of spws, because these have much noise. Only do this is
#       the previous dirty image is not just noise.

uvcontsub(vis='as2uds10.ms.split',              # Original file
          fitspw='2~13:5~58',                   # Deselect noise at spw edges
          combine='spw',                        # As we have left out entire spw's, set combine
          excludechans=True,                    # We want to exclude above chans
          want_cont=False)                      # Makes a continuum image

# Step 8: Load in physical parameters
# Expl: For the next step and the steps hereafter, we need the redshift and FWHM.
#       Find these using Birkin 2021 et al. table 2, open with file Read_Table_2_Birkin_2021


# Step 9: Cleaning
# Expl: Using the tclean method in CASA, we can try to suppress the noise around
#       the oject. We draw a mask of 6" using the interactive cleaning method and 
#       clean to a leven of 1sigma, because we have a low SNR.
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

# Step 10: 0-th moment map
# Expl: We make an integrated intensity map by collapsing the cube. This should
#       look something like Marta's right figures in her figure 2 of Frias Castillo
#       et al. 2023. For the channels to collapse we look at the FWHM found in step
#       6. Look at the velocities of the channels, we collapse on both sides for
#       the value of the FWHM (so 2sigma essentially).

immoments(imagename='as2uds10_128_05_100.split.cube.image',     # Original file
          chans='15~25',                        # Only collapse channels that fall in 2sigma FWHM
          moments=0,                            # Moment 0
          outfile='as2uds10_128_05_100.split.cube.image.mom0')  # Destination file

# Step 11: fits file of .mom0 map
# Expl: Make a file to read to make a figure of the .mom0 map

exportfits(imagename='as2uds10_128_05_100.split.cube.image.mom0',   # Original file
           fitsimage='as2uds10_128_05_100.split.cube.image.mom0.fits')  # Destination file

# Step 12: fits file of .cube map
# Expl: Make a file to read the RMS per channel using the .cube image

exportfits(imagename='as2uds10_128_05_100.split.cube.image',   # Original file
           fitsimage='as2uds10_128_05_100.split.cube.image.fits')  # Destination file

# Step 13: Get out the profile
# Expl: Imview the .mom0 file in imview spectral profile tool and load in the region
#       with 2.5 arcsec radius (already made). Center it on brightest the pixel
#       Set x to velocity and y to flux density. Then, save profile as .txt file. 
#       No function needed for this. 












