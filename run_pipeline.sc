#!/bin/bash -l

# Add CASA to path
export PATH=$PATH:'/home/sw/casa-5.4.1/bin/'

# Set MIRIAD environment
export MIR=/home/soft/miriad/
export MIRARCH=linux64
export MIRARCHD=$MIR/$MIRARCH
export MIRBIN=$MIRARCHD/bin
export MIRLIB=$MIRARCHD/lib
export MIRMAN=$MIRARCHD/man

export PATH=$MIRBIN:$PATH
export LD_LIBRARY_PATH=$MIRLIB:$LD_LIBRARY_PATH
export MANPATH=$MIRMAN:$MANPATH

export MIRCAT=$MIR/cat
export MIRPDOC=$MIR/doc
export MIRINC=$MIR/inc
export MIRPROG=$MIR/prog
export MIRSUBS=$MIR/subs

export MIRDEF=.
export PGPLOT_DIR=$MIRLIB
export PGPLOT_FONT=$MIRLIB/grfont.dat
export PGPLOT_RGB=$MIRLIB/rgb.txt

# Set path to mirpy
export MIRPY_DIR='/lemonpi/jallison/software/mirpy-0.2.0/'

# Add this directory to path
export PIPELINE='/lemonpi/jallison/thunderkat_hi_pipeline/'
export PATH=$PATH:$PIPELINE

# Set input file
input_file=$1

# Copy data to local directory
# casa -c $PIPELINE/split.py $input_file --nologger --log2term --nogui

# Add missing rest frequency field in MS table
# python2.7 $PIPELINE/add_restfreq.py $input_file

# Unflag data
# casa -c $PIPELINE/unflag.py $input_file --nologger --log2term --nogui

# Convert into UVFITS
# casa -c $PIPELINE/fits.py $input_file --nologger --log2term --nogui

# Run calibration and imaging
python2.7 $PIPELINE/cal_image.py --input_file $input_file

# Stack target and reference spectra
# python2.7 $PIPELINE/stack_spectra.py --input_file $input_file
