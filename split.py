# Import modules
import os
import sys
import numpy as np
import glob
from astropy.io import ascii

# Define and read input file
input_file = os.getcwd()+'/'+sys.argv[3]
input_lines = ascii.read(input_file, format='commented_header', comment='#', delimiter=' ')

# Loop over input lines
for line in input_lines:

	spw = ''
	if 'TRUE' in line['do_spw']:

		if line['res'] == '4k':

			# Define spectral properties of the data
			chunk_width = 20. # MHz
			chan_width = 208.984e-3 # MHz
			start_freq = 856. # MHz
			rest_freq = 1420.40575177 # MHz

		elif line['res'] == '32k':

			# Define spectral properties of the data
			chunk_width = 20. # MHz
			chan_width = 26.123046875e-3 # MHz
			start_freq = 856. # MHz
			rest_freq = 1420.40575177 # MHz

		# Calculate spectral window
		nchans = np.ceil(chunk_width/chan_width)
		chan_centre = np.round((rest_freq-start_freq)/chan_width)
		chan_low = np.floor(chan_centre - 0.5*(nchans))
		chan_high = np.ceil(chan_centre + 0.5*(nchans))
		spw =  '0:%d~%d'%(chan_low,chan_high)

	# Define paths to data
	in_ms = line['in_ms']
	out_path = os.getcwd()

	# Define observation ID and fields list
	obsid = line['id']
	fields = [line['targ_name'],line['pcal_name'],line['scal_name']]

	# Loop over objects and split
	for i in range(0,len(fields)):

		out_ms = out_path+'/%s/%s.ms'%(obsid,fields[i])

		# Define name of output MS and check for existence
		out_ms = out_path + '/%s/%s.ms'%(obsid,fields[i])
		if not os.path.exists('/'.join(out_ms.split('/')[0:-1])):
			os.system('mkdir -p %s'%'/'.join(out_ms.split('/')[0:-1]))

		if not os.path.exists(out_ms):
			mstransform(vis=in_ms,
				outputvis=out_ms,
				field=fields[i],
				spw=spw,
				datacolumn='data',
				timeaverage=False,
				chanaverage=False,
				usewtspectrum=True,
				realmodelcol=True)



