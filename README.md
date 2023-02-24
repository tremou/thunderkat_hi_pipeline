# thunderkat_hi_pipeline

The two txt files are for the source’s coords (source.txt) and the other is for the data information (j1848.txt) , resolution refant etc.. 
Both need to be edited in the beginning. 

The pipeline uses CASA 5.7.0-134 and python2.7 makes the following major steps: 

```source 01.run_pipeline_idia.sc j1848.txt```
- Copying part of the data locally (split) 
- Adding the rest frequencies in the tables 
- Unflagging any previous applied flags
- Converting the ms into fits that can be read by miriad. 

```source 02.run_pipeline_idia.sc j1848.txt```
- Bandpass calibration and imaging of the target 

```source 03.run_pipeline_idia.sc j1848.txt```
- Stacking the spectra (spectra.dat) and plotting them (pgplot) . It will assume by default the list of the sources in the source.txt 
