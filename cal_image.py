#!/usr/bin/env python2.7

## -- Initialize -- ##

# Import modules
import sys
import os
import numpy as np
sys.path.append(os.environ['MIRPY_DIR'])
from mirpy import miriad
from astropy.io import fits
from astropy.io import ascii
from astropy import wcs
from astropy import units as u
from astropy.coordinates import SkyCoord
import argparse

# Set user options
parser = argparse.ArgumentParser()
parser.add_argument('--input_file',type=str,help='set path to input file')
parser.add_argument('--source_list',default='./source_list.txt',type=str,help='set path to source list')
options = parser.parse_args()

# Define and read input file
input_file = os.getcwd()+'/'+options.input_file
input_lines = ascii.read(input_file, format='commented_header', comment='#', delimiter=' ')

# Define path to data
in_path = './' # os.getcwd()

# Define some constants
HI_FREQ = 1420.40575177 # MHz
LIGHT_SPEED = 2.99792458e5 # km/s

# Define whether fixed velocity bins
fixed_vel_bins = True

# Define whether HI emission filter
hi_filter = True

# Define method to find peak continuum in a small region
def find_peak(source_coord, image_name, box_size=11):

    # Open Image
    hdulist = fits.open(image_name)

    # Read data and header
    scidata = hdulist[0].data
    sci=np.squeeze(scidata)
    prihdr=hdulist[0].header

    # Make the header bi-dimensional so wcs does not get confused
    del prihdr['CTYPE4']
    del prihdr['CDELT4']
    del prihdr['CRVAL4']
    del prihdr['CRPIX4']
    del prihdr['CTYPE3']
    del prihdr['CDELT3']
    del prihdr['CRVAL3']
    del prihdr['CRPIX3']
    del prihdr['NAXIS3']
    del prihdr['NAXIS4']
    del prihdr['NAXIS']
    prihdr['NAXIS']=2

    # Load the WCS coordinate system module
    w=wcs.WCS(prihdr)
    px,py=w.wcs_world2pix(source_coord.ra.deg,source_coord.dec.deg,1)

    # Find position of peak flux density from region around source
    box = sci[int(round(py,0))-int(0.5*box_size)-1:int(round(py,0))+int(0.5*box_size),
         int(round(px,0))-int(0.5*box_size)-1:int(round(px,0))+int(0.5*box_size)]
    peak_flux = np.max(box)
    peak_pos = np.where(box == peak_flux)
    peak_py = int(round(py,0)) - int(0.5*box_size) + peak_pos[0]
    peak_px = int(round(px,0)) - int(0.5*box_size) + peak_pos[1]

    # Convert peak position to world coordinates
    peak_ra,peak_dec=w.wcs_pix2world(peak_px,peak_py,1)
    peak_coord = SkyCoord(ra=peak_ra,dec=peak_dec, unit=(u.degree, u.degree), frame='icrs')

    # Close file
    hdulist.close()

    # Return peak flux and position
    # print peak_flux, peak_px, peak_py
    return peak_flux, peak_coord

## -- Loop over input lines -- ##
for line in input_lines:

	# Define observation ID and objects list
	obsid = str(line['id'])
	names = [line['targ_name'],line['pcal_name'],line['scal_name']]

	# Define calibrators and target data
	bpass_cal = line['pcal_name']
	flux_cal = line['pcal_name']
	gain_cal = line['scal_name']
	target = line['targ_name']

	# Define channel range in which Galactic HI is likely to be located
	hi_chans = [int(line['hi_chans'].split(',')[0]),int(line['hi_chans'].split(',')[1])]

	# Define reference antenna
	refant = int(line['refant'])

	'''

	## -- Common data preparation -- ##
	for name in names:

		# Convert to MIRIAD UVFITS format
		in_file = in_path+'/'+obsid+'/'+name+'.fits'
		out_file = in_path+'/'+obsid+'/'+name+'.uv'
		if not os.path.exists(out_file):
			fits_input = {'In':in_file,
				'op':'uvin',
				'out':out_file,
				'velocity':'lsr'}
			miriad.fits(**fits_input)

		# Flag auto correlations
		in_file = in_path+'/'+obsid+'/'+name+'.uv'
		uvflag_input = {'vis':in_file,
			'select':'auto',
			'flagval':'flag'}
		miriad.uvflag(**uvflag_input)

		# Additional flags
		if os.path.exists(in_path+'/'+obsid+'/uvflag.txt'):
			flags = ascii.read(in_path+'/'+obsid+'/uvflag.txt', format='commented_header', comment='#')
			for flag in flags:
				uvflag_input = {'vis':in_file,
					'select':str(flag['select']),
					'flagval':'flag'}
				miriad.uvflag(**uvflag_input)

		# Add 1934-638 to header information for flux model
		if '193' in name:
			miriad.puthd(In=in_path+'/'+obsid+'/'+name+'.uv/source',
					value='%s'%('1934-638'))

		# Add HI rest frequency to header information
		miriad.puthd(In=in_path+'/'+obsid+'/'+name+'.uv/restfreq',
				value=HI_FREQ*1.e-3)


	## -- Calibrate bandpass calibrator data -- ##

	# Flag Galactic HI line channels
	uvflag_input = {'vis':in_path+'/'+obsid+'/'+bpass_cal+'.uv',
		'line':'channel,%d,%d,1,1'%(hi_chans[1]-hi_chans[0],hi_chans[0]),
		'flagval':'flag'}
	miriad.uvflag(**uvflag_input)

	# Use MFCAL to solve for bandpass and gains
	mfcal_input = {'vis':in_path+'/'+obsid+'/'+bpass_cal+'.uv',
		'stokes':'xx,yy',
		'refant':refant,
		'interval':'1,1,1e6',
		'options':'interpolate'}
	miriad.mfcal(**mfcal_input)

	# Smooth bandpass solution
	# gpedit_input = {'vis':in_path+'/'+obsid+'/'+bpass_cal+'.uv',
	#	'options':'hanning',
	#	'width':10}
	# miriad.gpedit(**gpedit_input)

	## -- Calibrate flux calibrator data -- ##
	if flux_cal != bpass_cal:

		# Flag Galactic HI line channels
		uvflag_input = {'vis':in_path+'/'+obsid+'/'+flux_cal+'.uv',
			'line':'channel,%d,%d,1,1'%(hi_chans[1]-hi_chans[0],hi_chans[0]),
			'flagval':'flag'}
		miriad.uvflag(**uvflag_input)

		# Use MFCAL to solve for bandpass and gains (note we don't use bandpass calibrator here)
		mfcal_input = {'vis':in_path+'/'+obsid+'/'+flux_cal+'.uv',
			'stokes':'xx,yy',
			'refant':refant,
			'interval':'1,1,1e6',
			'options':'interpolate'}
		miriad.mfcal(**mfcal_input)

		# Smooth bandpass solution
		# gpedit_input = {'vis':in_path+'/'+obsid+'/'+flux_cal+'.uv',
        #		'options':'hanning',
        #		'width':10}
		# miriad.gpedit(**gpedit_input)


	## -- Calibrate gain calibrator data -- ##
	if (gain_cal != flux_cal):

		if (gain_cal != bpass_cal):

			# Copy bandpass solution from bandpass calibrator
			gpcopy_input = {'vis':in_path+'/'+obsid+'/'+bpass_cal+'.uv',
				'out':in_path+'/'+obsid+'/'+gain_cal+'.uv',
				'mode':'copy',
				'options':'nocal,nopol'}
			miriad.gpcopy(**gpcopy_input)

			# Use MFCAL to solve for time varying gains 
			mfcal_input = {'vis':in_path+'/'+obsid+'/'+gain_cal+'.uv',
				'stokes':'xx,yy',
				'refant':refant,
				'interval':'5,5',
				'options':'nopassol'}
			miriad.mfcal(**mfcal_input)

		# Correct flux scale
		gpboot_input = {'vis':in_path+'/'+obsid+'/'+gain_cal+'.uv',
			'cal':in_path+'/'+obsid+'/'+flux_cal+'.uv'}
		miriad.gpboot(**gpboot_input)

		# Correct spectral slope based on flux cal model (probably overkill)
		# mfboot_input = {'vis':in_path+'/'+obsid+'/'+flux_cal+'.uv,'+in_path+'/'+obsid+'/'+gain_cal+'.uv',
		#	'select':'source(1934-638)'}
		# miriad.mfboot(**mfboot_input)


	## -- Copy and apply bandpass and time-varying gain solutions to the target data -- ##

	# Copy solutions from gain calibrator
	gpcopy_input = {'vis':in_path+'/'+obsid+'/'+gain_cal+'.uv',
				'out':in_path+'/'+obsid+'/'+target+'.uv',
				'mode':'copy'}
	miriad.gpcopy(**gpcopy_input)

	# Apply solutions
	in_file = in_path+'/'+obsid+'/'+target+'.uv'
	out_file = in_path+'/'+obsid+'/'+target+'.uv.cal'
	if not os.path.exists(out_file):
		uvcat_input = {'vis':infile,
				'out':out_file}
		miriad.uvcat(**uvcat_input)


	## -- Selfcal loop -- ##
	selfcal_intervals = [] # [2,2]
	selfcal_nfbins = [] # [1,1]
	selfcal_options = [] # ['phase','phase']
	selfcal_ind = 0
	for selfcal_ind in range(0,len(selfcal_intervals)):

		# image
		vis_file = in_path+'/'+obsid+'/'+target+'.uv.cal'
		map_file = in_path+'/'+obsid+'/'+target+'.mfs.imap.cal%d' % (selfcal_ind)
		beam_file = in_path+'/'+obsid+'/'+target+'.mfs.ibeam.cal%d' % (selfcal_ind)
		if not os.path.exists(map_file):
			invert_input = {'vis':vis_file,
				'map':map_file,
				'beam':beam_file,
				'imsize':'4096,4096',
				'cell':'4,4,res',
				'robust':'-0.5',
				'stokes':'ii',
				'options':'double,mfs',
				'mode':'fft',
				'slop':'1,interpolate'}
			miriad.invert(**invert_input)

		# calculate image noise
		sigest_input = {'In':map_file,
			'region':'box(0,0,128,128)'}
		sigest_output = miriad.sigest(**sigest_input)
		image_noise = float(sigest_output.split('\n')[-1].split(' ')[-1])

		# clean
		out_file = in_path+'/'+obsid+'/'+target+'.mfs.icmp.cal%d' % (selfcal_ind)
		if not os.path.exists(out_file):
			clean_input = {'map':map_file,
				'beam':beam_file,
				'out':out_file,
				'gain':'0.01',
				'options':'negstop',
				'cutoff':'%f'%(3*image_noise),
				'niters':'1e6',
				'speed':'0'}
			miriad.clean(**clean_input)

		# restor
		out_file = in_path+'/'+obsid+'/'+target+'.mfs.icln.cal%d' % (selfcal_ind)
		model_file = in_path+'/'+obsid+'/'+target+'.mfs.icmp.cal%d' % (selfcal_ind)
		if not os.path.exists(out_file):
			restor_input = {'model':model_file,
				'beam':beam_file,
				'map':map_file,
				'out':out_file}
			miriad.restor(**restor_input)

		# selfcal
		vis_file = in_path+'/'+obsid+'/'+target+'.uv.cal'
		selfcal_input = {'vis':vis_file,
			'model':model_file,
			'interval':'%f'%(selfcal_intervals[selfcal_ind]),
			'nfbin':'%f'%(selfcal_nfbins[selfcal_ind]),
			'options':'%s'%(selfcal_options[selfcal_ind])}
		miriad.selfcal(**selfcal_input)

		if (selfcal_ind == len(selfcal_intervals)-1):
			
			# Apply selfcal solutions
			vis_file = in_path+'/'+obsid+'/'+target+'.uv.cal'
			out_file = in_path+'/'+obsid+'/'+target+'.uv.cal.tmp'
			if not os.path.exists(out_file):
				uvcat_input = {'vis':vis_file,
					'out':out_file}
				miriad.uvcat(**uvcat_input)
				os.system('rm -rf %s' % (vis_file))
				os.system('rsync -P -rte ssh %s/* %s' % (out_file, vis_file))
				os.system('rm -rf %s' % (out_file))


	## -- Final continuum image -- ##

	# image
	vis_file = in_path+'/'+obsid+'/'+target+'.uv.cal' 
	map_file = in_path+'/'+obsid+'/'+target+'.mfs.imap.cal%d' % (selfcal_ind+1)
	beam_file = in_path+'/'+obsid+'/'+target+'.mfs.ibeam.cal%d' % (selfcal_ind+1)
	if not os.path.exists(map_file):
		invert_input = {'vis':vis_file,
					'map':map_file,
					'beam':beam_file,
					'imsize':'4096,4096',
					'cell':'4,4,res',
					'robust':'-0.5',
					'stokes':'ii',
					'options':'double,mfs',
					'mode':'fft',
					'slop':'1,interpolate'}
		miriad.invert(**invert_input)

	# calculate image noise
	sigest_input = {'In':map_file,
		'region':'box(0,0,128,128)'}
	sigest_output = miriad.sigest(**sigest_input)
	image_noise = float(sigest_output.split('\n')[-1].split(' ')[-1])

	# clean
	out_file = in_path+'/'+obsid+'/'+target+'.mfs.icmp.cal%d' % (selfcal_ind+1)
	if not os.path.exists(out_file):
		clean_input = {'map':map_file,
				'beam':beam_file,
				'out':out_file,
				'gain':'0.01',
				'options':'negstop',
				'cutoff':'%f'%(3.*image_noise),
				'niters':'1e6',
				'speed':'0'}
		miriad.clean(**clean_input)

	# restor
	out_file = in_path+'/'+obsid+'/'+target+'.mfs.icln.cal%d' % (selfcal_ind+1)
	model_file = in_path+'/'+obsid+'/'+target+'.mfs.icmp.cal%d' % (selfcal_ind+1)
	if not os.path.exists(out_file):
		restor_input = {'model':model_file,
				'beam':beam_file,
				'map':map_file,
				'out':out_file}
		miriad.restor(**restor_input)

	'''

	## -- Loop over source list  --- ##

	# Since MIRIAD does not correct for the w-term when imaging, we
	# phase rotate the visibilities to the position of each source and
	# produce an image cube and spectrum

	# Read source list
	sources = ascii.read(options.source_list, format='commented_header', comment='#', delimiter=' ')

	source_index = 0
	for source in sources:

		# Create source position object
		source_coord = SkyCoord(ra=source['ra'],dec=source['dec'],unit=(u.hourangle, u.deg),frame='icrs')

		# Create source name
		if source_coord.dec.degree < 0.:
        	        source_dec_sign = '-'
	        else:
	                source_dec_sign = '+'
		source_name = 'J%02d%02d%05.2f%s%02d%02d%04.1f' % (int(np.abs(source_coord.ra.hms[0])),
                                int(np.abs(source_coord.ra.hms[1])),
                                np.abs(source_coord.ra.hms[2]),
				str(source_dec_sign),
                                int(np.abs(source_coord.dec.dms[0])),
                                int(np.abs(source_coord.dec.dms[1])),
                                np.abs(source_coord.dec.dms[2]))

		# phase rotate calibrated visibilities to source position
		vis_file = in_path+'/'+obsid+'/'+target+'.uv.cal'
		out_file = in_path+'/'+obsid+'/'+source_name+'.uv'
		if not os.path.exists(out_file):
			uvedit_ra = '%02d,%02d,%f' % (int(np.abs(source_coord.ra.hms[0])),
                                int(np.abs(source_coord.ra.hms[1])),
                                np.abs(source_coord.ra.hms[2]))
			uvedit_dec = '%s%02d,%02d,%f' % (str(source_dec_sign),
                                int(np.abs(source_coord.dec.dms[0])),
                                int(np.abs(source_coord.dec.dms[1])),
                                np.abs(source_coord.dec.dms[2]))
			uvedit_input = {'vis':vis_file,
					'ra':uvedit_ra,
					'dec':uvedit_dec,
					'out':out_file}
			miriad.uvedit(**uvedit_input)

		## -- Continuum Imaging of source  --- ##

		# continuum image
		vis_file = in_path+'/'+obsid+'/'+source_name+'.uv' 
		map_file = in_path+'/'+obsid+'/'+source_name+'.mfs.imap'
		beam_file = in_path+'/'+obsid+'/'+source_name+'.mfs.ibeam'
		if not os.path.exists(map_file):
			invert_input = {'vis':vis_file,
					'map':map_file,
					'beam':beam_file,
					'imsize':'4096,4096',
					'cell':'4,4,res',
					'robust':'0.5',
					'stokes':'ii',
					'options':'double,mfs',
					'mode':'fft',
					'slop':'1,interpolate'}
			miriad.invert(**invert_input)

		# calculate image noise
		sigest_input = {'In':map_file,
			'region':'box(0,0,128,128)'}
		sigest_output = miriad.sigest(**sigest_input)
		image_noise = float(sigest_output.split('\n')[-1].split(' ')[-1])

		# clean
		out_file = in_path+'/'+obsid+'/'+source_name+'.mfs.icmp'
		if not os.path.exists(out_file):
			clean_input = {'map':map_file,
					'beam':beam_file,
					'out':out_file,
					'gain':'0.01',
					'options':'negstop',
					'cutoff':'%f'%(3.*image_noise),
					'niters':'1e6',
					'speed':'0'}
			miriad.clean(**clean_input)

		# restor
		out_file = in_path+'/'+obsid+'/'+source_name+'.mfs.icln'
		model_file = in_path+'/'+obsid+'/'+source_name+'.mfs.icmp'
		if not os.path.exists(out_file):
			restor_input = {'model':model_file,
					'beam':beam_file,
					'map':map_file,
					'out':out_file}
			miriad.restor(**restor_input)

		# convert to fits
		in_file = in_path+'/'+obsid+'/'+source_name+'.mfs.icln'
		out_file = in_path+'/'+obsid+'/'+source_name+'.mfs.icln.fits'
		if not os.path.exists(out_file):
			fits_input = {'In':in_file,
					'op':'xyout',
					'out':out_file}
			miriad.fits(**fits_input)

		# Find position of source peak continuum
		in_file = in_path+'/'+obsid+'/'+source_name+'.mfs.icln.fits'
		peak_flux, peak_coord = find_peak(source_coord,in_file,box_size=11)


		## -- Continuum Subtraction -- ##

		# uvmodel
		vis_file = in_path+'/'+obsid+'/'+source_name+'.uv'
		model_file = in_path+'/'+obsid+'/'+source_name+'.mfs.icmp'
		out_file = in_path+'/'+obsid+'/'+source_name+'.uv.uvmodel'
		if not os.path.exists(out_file):
			uvmodel_input = {'vis':vis_file,
					'model':model_file,
					'options':'subtract,mfs',
					'out':out_file}
			if fixed_vel_bins:
					uvmodel_input['line'] = 'velocity,75,-200,5.6'
			miriad.uvmodel(**uvmodel_input)


		# uvlin
		vis_file = in_path+'/'+obsid+'/'+source_name+'.uv.uvmodel'
		out_file = in_path+'/'+obsid+'/'+source_name+'.uv.uvlin'
		if not os.path.exists(out_file):
			uvlin_input = {'vis':vis_file,
					'out':out_file,
					'order':2,
					'mode':'line',
					'options':'nowindow'}
			if fixed_vel_bins:
					uvlin_input['line'] = 'velocity,75,-200,5.6'
					uvlin_input['chans'] = '0,18,54,1e9'
			else:
					uvlin_input['chans'] = '0,%d,%d,1e9' % (hi_chans[0],hi_chans[1]),
			miriad.uvlin(**uvlin_input)



		## -- Spectral-line Imaging of source  --- ##		

		# invert
		vis_file = in_path+'/'+obsid+'/'+source_name+'.uv.uvlin'
		map_file = in_path+'/'+obsid+'/'+source_name+'.cube.imap'
		beam_file = in_path+'/'+obsid+'/'+source_name+'.cube.ibeam'
		if not os.path.exists(map_file):
			invert_input = {'vis':vis_file,
					'map':map_file,
					'beam':beam_file,
					'imsize':'2048,2048',
					'cell':'4,4,res',
					'robust':'0.5',
					'stokes':'ii',
					'options':'double',
					'mode':'fft',
					'slop':'1,zero'}
			if fixed_vel_bins:
					invert_input['line'] = 'velocity,75,-200,5.6'
			if hi_filter:
					invert_input['select'] = '-uvrange(0,5)'
			miriad.invert(**invert_input)

		# calculate image noise
		sigest_input = {'In':map_file,
		        'region':'box(0,0,128,128)'}
		sigest_output = miriad.sigest(**sigest_input)
		image_noise = float(sigest_output.split('\n')[-1].split(' ')[-1])

		# Make continuum mask for cleaning HI absorption
		in_file = in_path+'/'+obsid+'/'+source_name+'.mfs.icln'
		out_file = in_path+'/'+obsid+'/'+source_name+'.mfs.mask'
		if not os.path.exists(out_file):
			maths_input = {'exp':'<%s>'%(in_file),
						'mask':'<%s>.gt.(%.8f)'%(in_file, 10.*image_noise),
						'region':'perc(100)',
						'out':out_file}
			miriad.maths(**maths_input)

		# alter third axis of mask to span the observed bandwidth
		in_file = in_path+'/'+obsid+'/'+source_name+'.mfs.mask'
		puthd_input = {'in':'%s/cdelt3'%(in_file),
					'value':'999'}
		miriad.puthd(**puthd_input)

		# regrid continuum mask to dimensions of cube
		in_file = in_path+'/'+obsid+'/'+source_name+'.mfs.mask'
		tin_file = in_path+'/'+obsid+'/'+source_name+'.cube.imap'
		out_file = in_path+'/'+obsid+'/'+source_name+'.mfs.mask.regrid'
		mask_success = False
		if not os.path.exists(out_file):
			regrid_input = {'in':in_file,
							'out':out_file,
							'tin':tin_file}
			try:
				miriad.regrid(**regrid_input)
				mask_success = True
			except:
				# Masking hasn't worked, likely continuum lower than cube noise threshold
				mask_success = False


		# clean
		out_file = in_path+'/'+obsid+'/'+source_name+'.cube.icmp'
		mask_file = in_path+'/'+obsid+'/'+source_name+'.mfs.mask.regrid'
		if not os.path.exists(out_file):
			clean_input = {'map':map_file,
					'beam':beam_file,
					'out':out_file,
					'gain':'0.01',
					'cutoff':'%f'%(5.*image_noise),
					'niters':'1e6',
					'speed':'0'}
			if mask_succes:
				clean_input['region'] = 'mask(%s)'%(mask_file)
			miriad.clean(**clean_input)

		# restor
		out_file = in_path+'/'+obsid+'/'+source_name+'.cube.icln'
		model_file = in_path+'/'+obsid+'/'+source_name+'.cube.icmp'
		if not os.path.exists(out_file):
			restor_input = {'model':model_file,
					'beam':beam_file,
					'map':map_file,
					'out':out_file}
			miriad.restor(**restor_input)


		## -- Extract spectra -- ##

		# It is important to note that one should be cautious when extracting
		# spectral line data from MIRIAD cubes using non-MIRIAD tasks. This
		# is due to the pixel sizes not always being a fixed angular size

		# Extract spectrum using MBSPECT
		in_file = in_path+'/'+obsid+'/'+source_name+'.cube.icln'
		out_file = in_path+'/'+obsid+'/'+source_name+'.mbspect'
		if peak_coord.dec.degree < 0.:
			peak_dec_sign = '-'
		else:
			peak_dec_sign = '+'
		mbspect_ra = '%02d:%02d:%f'%(int(np.abs(peak_coord.ra.hms[0])),
					int(np.abs(peak_coord.ra.hms[1])),
					np.abs(peak_coord.ra.hms[2]))
		mbspect_dec = '%s%02d:%02d:%f'%(str(peak_dec_sign),
					int(np.abs(peak_coord.dec.dms[0])),
					int(np.abs(peak_coord.dec.dms[1])),
					np.abs(peak_coord.dec.dms[2]))
		if not os.path.exists(out_file):
			mbspect_input = {'In':in_file,
					'out':out_file,
					'coord':'%s,%s' % (mbspect_ra,mbspect_dec),
					'width':'1,1',
					'xaxis':'felo',
					'yaxis':'average'}
			miriad.mbspect(**mbspect_input)

		# Write spectrum to ascii file using IMSPECT
		in_file = in_path+'/'+obsid+'/'+source_name+'.mbspect'
		out_file = in_path+'/'+obsid+'/'+source_name+'.imspect.log'
		if not os.path.exists(out_file):
			imspect_input = {'In':in_file,
					'hann':1,
					'log':out_file}
			miriad.imspect(**imspect_input)

		# Read spectrum in astropy table
		imspect_log = in_path+'/'+obsid+'/'+source_name+'.imspect.log'
		spectrum = ascii.read(imspect_log,data_start=3,format='no_header')
		spectrum.rename_column('col1', 'index')
		spectrum.rename_column('col2', 'radio_vel[km/s]')
		spectrum.rename_column('col3', 'flux[Jy]')

		# Estimate noise spectrum from cube
		in_file = in_path+'/'+obsid+'/'+source_name+'.cube.icln'
		imstat_log = in_path+'/'+obsid+'/'+source_name+'.cube.imstat.log'
		imstat_input = {'In':in_file,
				'region':'box(0,0,128,128)',
				'plot':'rrms',
				'options':'guaranteespaces',
				'log':imstat_log}
		miriad.imstat(**imstat_input)
		imstat_data=ascii.read(imstat_log,header_start=7,data_start=8,data_end=-5)
		spectrum['flux_noise[Jy]'] = imstat_data['Rrms']

		# Calculate other columns for output spectrum
		spectrum['freq[MHz]'] = HI_FREQ*(1.-spectrum['radio_vel[km/s]']/LIGHT_SPEED)
		spectrum['cont[Jy]'] = peak_flux
		spectrum['abs'] = spectrum['flux[Jy]']/spectrum['cont[Jy]']
		spectrum['abs_noise'] = spectrum['flux_noise[Jy]']/spectrum['cont[Jy]']

		# Write out spectrum to file
		ascii.write(spectrum['freq[MHz]','radio_vel[km/s]','flux[Jy]',
        		                'flux_noise[Jy]','cont[Jy]',
	        	                'abs','abs_noise'],
			in_path+'/'+obsid+'/'+'%s_spectrum.dat' % (source_name),
			formats={'freq[MHz]':'%.8e','radio_vel[km/s]':'%.8e',
				'flux[Jy]':'%.8e','flux_noise[Jy]':'%.8e',
				'cont[Jy]':'%.8e','abs':'%.8e','abs_noise':'%.8e'},
			format='commented_header', comment='#',
				delimiter=' ',overwrite=True)

		# Increment source index
		source_index += 1
