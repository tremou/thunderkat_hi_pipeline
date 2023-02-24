# thunderkat_hi_pipeline

The two txt files are for the source’s coords (source list ) and the other is for the ms info , resolution refant etc.. 
Both need to be edited in the beginning. 

The pipeline uses CASA 5.7.0-134 and python2.7 makes the following major steps: 

>source 01.run_pipeline_idia.sc test.txt
- Copies part of the data locally (split) 
- Adding the rest frequencies in the tables 
- Unflags any previous flags
- Converting the ms into fits that can be read by miriad. 

>source 02.run_pipeline_idia.sc test.txt
- Bandpass calibration and image the target 

>source 03.run_pipeline_idia.sc source_list.txt
- Stack spectra and plots them (pgplot) 
