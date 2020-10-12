## Copyright James Allison, 2019, All rights reserved

## Stack spectra

## -- Initialize -- ##

# Import required modules
import os
import sys
import numpy as np
import glob
from astropy.io import ascii
from astropy.table import Table, vstack
import argparse

# Import matplotlib and set properties
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rc
import matplotlib.gridspec as gridspec
# rc('text', usetex=True)
# rc('font',**{'family':'serif','serif':['serif'],'size':20})
matplotlib.rcParams['xtick.direction'] = 'in'
matplotlib.rcParams['ytick.direction'] = 'in'

# Set user options
parser = argparse.ArgumentParser()
parser.add_argument('--input_file',type=str,help='set path to input file')
parser.add_argument('--source_list',default='./source_list.txt',type=str,help='set path to source list')
parser.add_argument('--vel_min',default=-790.,type=float,help='set minimum velocity for plotting')
parser.add_argument('--vel_max',default=790.,type=float,help='set maximum velocity for plotting')
options = parser.parse_args()

# Define and read input file
input_file = os.getcwd()+'/'+options.input_file
input_lines = ascii.read(input_file, format='commented_header', comment='#', delimiter=' ')

# Define and read source list
source_list = os.getcwd()+'/'+options.source_list
sources = ascii.read(source_list, format='commented_header', comment='#', delimiter=' ')
sources['name'] = None
ind = 0
for source in sources:
	sources['name'][ind] = 'J'+source['ra'].replace(':','')[0:6]+source['dec'].replace(':','')[0:7]
	ind += 1

# Define constants
LIGHT_SPEED = 2.99792458e5 # km/s
HI_FREQ = 1420.40575177 # MHz

# List all spectra
spectra = glob.glob("%s/*/*_spectrum.dat" % os.getcwd())

# Loop over spectra and construct composite target and reference spectrum
target_spectrum = Table(names=('radio_vel[km/s]','abs','abs_noise'))
reference_spectrum = Table(names=('radio_vel[km/s]','abs','abs_noise'))
for spectrum in spectra:

	spectrum_name = spectrum.split('/')[-1].split('_spectrum.dat')[0]
	if '+' in spectrum_name:
		spectrum_sign = '+'
	else:
		spectrum_sign = '-'
	spectrum_ra = spectrum_name.split('J')[-1].split(spectrum_sign)[0]
	spectrum_dec = spectrum_name.split('J')[-1].split(spectrum_sign)[1]
	spectrum_name = 'J'+spectrum_ra[0:6]+spectrum_sign+spectrum_dec[0:6]

	reference = sources['reference'][sources['name'] == spectrum_name]

	if len(reference):

		if reference[0] == 'FALSE':

			target_spectrum =  vstack([target_spectrum,ascii.read(spectrum)])

			target_name = spectrum_name

		else:

			reference_spectrum = vstack([reference_spectrum,ascii.read(spectrum)])

if len(reference_spectrum) == 0:

	reference_spectrum = target_spectrum

## -- Bin composite spectra -- ##

# Construct velocity bins
dvel = np.abs(np.median(np.diff(target_spectrum['radio_vel[km/s]'])))
bin_edges = np.arange(np.min(target_spectrum['radio_vel[km/s]'])-dvel*0.5, np.max(target_spectrum['radio_vel[km/s]'])+0.5*dvel, dvel)


null_array = np.nan*np.zeros(len(bin_edges))
target_spectrum_binned = Table([null_array,null_array,null_array],names=('radio_vel[km/s]','abs','abs_noise'))
reference_spectrum_binned = Table([null_array,null_array,null_array],names=('radio_vel[km/s]','abs','abs_noise'))
ind = 0
for bin_edge in bin_edges-1:

	# Target spectrum

	# Set up bin filter
	truth = (target_spectrum['radio_vel[km/s]']>=bin_edge) & (target_spectrum['radio_vel[km/s]']<bin_edge+dvel) \
			& (target_spectrum['abs_noise']>0.0) & (np.isfinite(target_spectrum['abs_noise'])) \
			& (~np.isnan(target_spectrum['abs_noise']))
	if np.size(truth[truth]) == 0:
		ind += 1
		continue

	# Calculate mean bin velocity, optical depth and noise
	weights = np.power(target_spectrum['abs_noise'][truth], -2.)
	sum_weights = np.sum(weights)
	target_spectrum_binned['radio_vel[km/s]'][ind] = np.sum(target_spectrum['radio_vel[km/s]'][truth]*weights)/sum_weights
	target_spectrum_binned['abs'][ind] = np.sum(target_spectrum['abs'][truth]*weights)/sum_weights
	target_spectrum_binned['abs_noise'][ind] = np.sqrt(np.sum(np.power(target_spectrum['abs_noise'][truth],2)*
								np.power(weights, 2)))/sum_weights

	# Reference spectrum

	# Set up bin filter
	truth = (reference_spectrum['radio_vel[km/s]']>=bin_edge) & (reference_spectrum['radio_vel[km/s]']<bin_edge+dvel) \
			& (reference_spectrum['abs_noise']>0.0) & (np.isfinite(reference_spectrum['abs_noise'])) \
			& (~np.isnan(reference_spectrum['abs_noise']))
	if np.size(truth[truth]) == 0:
		ind += 1
		continue

	# Calculate mean bin velocity, optical depth and noise
	weights = np.power(reference_spectrum['abs_noise'][truth], -2.)
	sum_weights = np.sum(weights)
	reference_spectrum_binned['radio_vel[km/s]'][ind] = np.sum(reference_spectrum['radio_vel[km/s]'][truth]*weights)/sum_weights
	reference_spectrum_binned['abs'][ind] = np.sum(reference_spectrum['abs'][truth]*weights)/sum_weights
	reference_spectrum_binned['abs_noise'][ind] = np.sqrt(np.sum(np.power(reference_spectrum['abs_noise'][truth],2)*
								np.power(weights, 2)))/sum_weights

	# Increment index
	ind += 1

## -- Write stacked spectra to file -- ##
ascii.write(target_spectrum_binned, './target_spectrum_binned.dat', 
			formats={'radio_vel[km/s]':'%.8e', 'abs':'%.8e','abs_noise':'%.8e'}, 
			delimiter=' ', format='commented_header', comment='#', overwrite=True)
ascii.write(reference_spectrum_binned, './reference_spectrum_binned.dat', delimiter=' ',  
			formats={'radio_vel[km/s]':'%.8e', 'abs':'%.8e','abs_noise':'%.8e'}, 
			format='commented_header', comment='#', overwrite=True)

## -- Make plot of binned spectra - ##
plt.ioff()
fig = plt.figure(figsize = (8,6))
#plt.tight_layout()
nxplots = 1
nyplots = 1
gs = gridspec.GridSpec(nyplots,nxplots)
gs.update(wspace=0.0, hspace=0.0)

# Initialize subplot
ax1 = plt.subplot(gs[0])

# Set limits
x_min = options.vel_min
x_max = options.vel_max
y_min = 1.5*100.*np.nanmin(reference_spectrum_binned['abs'])
y_max = 1.5*100.*np.nanmax(reference_spectrum_binned['abs'])
# y_max = 3.*100.*np.nanmax(target_spectrum_binned['abs'])
# offset = 3.*100.*(np.nanmin(target_spectrum_binned['abs'])-np.nanmax(reference_spectrum_binned['abs']))
offset = 0.
y_min += offset
ax1.set_xlim(x_min, x_max)
ax1.set_ylim(y_min, y_max)

# Do plotting
h1 = ax1.step(target_spectrum_binned['radio_vel[km/s]'], 100.*target_spectrum_binned['abs'], where='mid', linestyle='-',
				color='r',linewidth=2., zorder=2, label=r'%s'%target_name)
ax1.fill_between(target_spectrum_binned['radio_vel[km/s]'], 0.0, 3.*100.*target_spectrum_binned['abs_noise'],
				facecolor='r', edgecolor='none', zorder=0, alpha=0.25)
ax1.fill_between(target_spectrum_binned['radio_vel[km/s]'], 0.0, -3.*100.*target_spectrum_binned['abs_noise'],
				facecolor='r', edgecolor='none', zorder=0, alpha=0.25)
h2 = ax1.step(reference_spectrum_binned['radio_vel[km/s]'], 100.*reference_spectrum_binned['abs']+offset, where='mid', linestyle='-',
 				color='b',linewidth=2., zorder=2, label=r'Reference')
ax1.fill_between(reference_spectrum_binned['radio_vel[km/s]'], offset, 3.*100.*reference_spectrum_binned['abs_noise']+offset,
				facecolor='b', edgecolor='none', zorder=0, alpha=0.25)
ax1.fill_between(reference_spectrum_binned['radio_vel[km/s]'], offset, -3.*100.*reference_spectrum_binned['abs_noise']+offset,
				facecolor='b', edgecolor='none', zorder=0, alpha=0.25)
ax1.axhline(0., color='k', linestyle=':', linewidth=1.0, zorder=1)
ax1.axhline(offset, color='k', linestyle=':', linewidth=1.0, zorder=1)
ax1.axvline(0., color='k', linestyle=':', linewidth=1.0, zorder=1)

# Set labels
ax1.set_xlabel(r"LSR velocity [km/s]", fontsize = 20)
ax1.set_ylabel(r"Fractional absorption [per cent]", fontsize = 20)

# Add legend
hs = h1+h2
labs = [h.get_label() for h in hs]
ax1.legend(loc='lower left',frameon=False,fontsize=16)

# Set minor ticks
ax1.minorticks_on()

# Set tick lengths
ax1.tick_params(bottom=True, left=True, top=True, right=True,
				length=10, width=1, which='major', direction='in')
ax1.tick_params(bottom=True, left=True, top=True, right=True,
				length=5, width=1, which='minor', direction='in')

# Write spectrum to PDF file
plt.savefig('./spectrum_binned.pdf')


