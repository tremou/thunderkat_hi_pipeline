# Import modules
import sys
from astropy.io import ascii

# Define and read input file
input_file = os.getcwd()+'/'+sys.argv[3]
input_lines = ascii.read(input_file, format='commented_header', comment='#', delimiter=' ')

# Define path to data
in_path = os.getcwd()

# Loop over input lines
for line in input_lines:

	# Define observation ID and fields list
	obsid = line['id']
	fields = [line['targ_name'],line['pcal_name'],line['scal_name']]

	# Loop over targes and split
	for i in range(0,len(fields)):

		# Define input MS
		in_ms = in_path+'/%s/%s.ms' % (obsid, fields[i])

		# Convert to fits
		fits = in_ms.replace('.ms','.fits')
		exportuvfits(vis=in_ms,
					fitsfile=fits,
					datacolumn='data')


