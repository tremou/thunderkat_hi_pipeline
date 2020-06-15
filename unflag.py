# Import modules
import sys
import os
from astropy.io import ascii

# Define and read input file
input_file = sys.argv[3]
input_lines = ascii.read(input_file, format='commented_header', comment='#', delimiter=' ')

# Define path to data
in_path = os.getcwd()

# Loop over input lines and unflag data
for line in input_lines:

	# Check unflag option
	do_unflag=line['do_unflag']
	if 'FALSE' in do_unflag:
		continue

	# Define observation ID and fields list
	obsid = line['id']
	fields = [line['targ_name'],line['pcal_name'],line['scal_name']]

	# Loop over targes and split
	for i in range(0,len(fields)):

		# Define input MS
		in_ms = in_path+'/%s/%s.ms' % (obsid, fields[i])

		# Unflag data
		flagdata(vis = in_ms, mode = 'unflag')
