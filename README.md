# thunderkat_hi_pipeline

The two txt files are for the source’s coords (source.txt) and the other (j1848.txt) are the inputs for the data information, resolution, refant etc.. 
Both need to be edited in the beginning accordingly.

The pipeline uses CASA version 5.7.0-134, python2.7 and MIRIAD. 
The three following major steps need to be run seperately. 

```source 01.run_pipeline_idia.sc j1848.txt```

It calls the split.py, add_restfreq.py, unflag.py, and fits.py python scripts performing the following steps: 
- Copying part of the data locally (split) 
- Adding the rest frequencies in the tables 
- Unflagging any previous applied flags
- Converting the ms into fits that can be read by miriad. 

The expected output products are measurement sets and uvfits files of each field.  

```source 02.run_pipeline_idia.sc j1848.txt```

It calls the cal_image.py python script performing the following steps:
- Data pre-processing. 
  - region definition to find he peak continuum position and to convert it into WCS coordinates 
  - defining the fields/names of the calibrators and target
  - setting reference antenna 
  - basic flagging
- Bandpass and gain calibration, and imaging (clean). 
  - bandpass calibration using the flux calibrator (primary calibrator) in order to correct as a function of frequency. 
  - solving for time varying gains using the secondary calibrator (phase calibrator). 
  - flux scaling
  - transfer and applying those corrections to the target field. 
  - phase rotation of the calibrated visibilities to the source position lsited in the source_list.txt file (MIRIAD does not account for w-term). 
  - continuum imaging (clean) and self calibration of the target and the list of teh sources listed in the source_list.txt file. 
  - continuum subtraction
  - spectral line imaging
  - spectra extraction into an ascii file (.dat)
  
The expected output products are freqeuncy cubes, images of the target and reference sources and the spectra ascii file. 


```source 03.run_pipeline_idia.sc j1848.txt```
- Stacking the spectra (*.dat) and plotting them (pgplot) . It will assume by default the list of the sources in the source_list.txt

The expected output products are the spectrum of the refernce source and the target (ascii files- reference_spectrum_binned.dat, target_spectrum_binned.dat) and a plot of the binned spectra (spectrum_binned.pdf). 
