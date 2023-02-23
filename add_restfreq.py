# Import modules
import os
import sys
#from pyrap.tables import table
from casacore.tables.table import table
from astropy.io import ascii

# Define and read input file
input_file = os.getcwd()+'/'+sys.argv[1]
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

		# Define rest frequency dictionary
		desc = {'name': 'REST_FREQUENCY',
			'_c_order': True,
			'comment': 'Line rest frequency',
			'dataManagerGroup': 'StandardStMan',
			'dataManagerType': 'StandardStMan',
			'keywords': {'MEASINFO': {'Ref': 'LSRK', 'type': 'frequency'},
			'QuantumUnits': ['Hz']},
			'maxlen': 0,
			'ndim': -1,
			'option': 0,
			'valueType': 'double'}

		# Add rest frequency field to table
		tt = table('%s/SOURCE'%(in_ms), readonly=False)
		if 'REST_FREQUENCY' not in tt.colnames():
			tt.addcols(desc)
		tt.done()
