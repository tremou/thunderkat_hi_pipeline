#!/bin/bash -l

# Set MIRIAD environment
export MIR="/scratch3/users/tremou/soft/miriad/"
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
export MIRPY_DIR='/scratch3/users/tremou/soft/mirpy-master/mirpy'

# Add this directory to path
export PIPELINE='/scratch3/users/tremou/soft/thunderkat_hi_pipeline-idia_j1848/'
export PATH=$PATH:$PIPELINE

# Set input file
input_file=$1

# Set up container aliases
alias casa_container="singularity exec /idia/software/containers/casa-stable.img"

# Copy data to local directory
casa_container casa -c $PIPELINE/split.py $input_file --nologger --log2term --nogui

# Add missing rest frequency field in MS table
casa_container python2.7 $PIPELINE/add_restfreq.py $input_file

# Unflag data
casa_container casa -c $PIPELINE/unflag.py $input_file --nologger --log2term --nogui

# Convert into UVFITS
casa_container casa -c $PIPELINE/fits.py $input_file --nologger --log2term --nogui

