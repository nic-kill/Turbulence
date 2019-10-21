"""
5/23/19
This script will perform velocity channel analysis (VCA) on an input PPV cube/FITS file. The spectral axis of the input cube will be down 
sampled with a Gaussian kernel in order to change the velocity resolution according to the user provided number of VCA bins (N_VCA). For example, 
if the input cube has a 100 channels --- each 1 km/s in width --- and the user requests 20 VCA bins, the spectral axis will be downsampled
such that the first returned VCA slope is the average of SPS slopes computed at every 1 km/s channel, the second VCA slope is the average SPS slope computed at ever channel
downsampled to 6.21 km/s resolution, the third ... at 11.42 km/s resolution, ... until the SPS is computed with the spectral axis downsampled to a single 100 km/s channel, 
which is essentially the integrated intensity image. The formula for computing the subsequent thickness of each velocity slice is: abs(vN - v0) / (N_VCA - 1), where vN is the velocity 
of the last spectral channel, v0 is the velocity of the first spectral channel. This is designed such that the last down-sample is simply the mean intensity along the spectral axis. 

The SPS slopes are computed utilizing a custom code called SPS_2D.py. The user specifies number of bins (spaced evenly in log space) used to compute the SPS of each velocity channel (N_bins).
See the header of that script for explicit details on the computation and uncerainties. A plot of the weighted average of the SPS slopes as a function of velocity thickness will be produced. 
A final input signifies whether the final plot and variables will be saved to disk in the form of a pdf and binary file, respectively. Note that the input cube must have beam and
spectral axis information defined in the header; additionally, the spectral axis is assumed to be in units of radial velocity [km/s]. 
A summary of the inputs:
Inputs:
1 - path to FITS image
2 - number of VCA bins
3 - number of annuli/bins used to compute SPS slope (see header of SPS_2D.py)
4 - 1 or 0 (True or False) to save output plot and variables in binary pickle file 
Usage:
ipython VCA.py /path/to/FITS/image/ 20 12 1
Example:
ipython SPS_2D.py SMC_ASKAPsoft_JD_2000min_5maj_Parkes_K_MOM0_ATCAREGRID.FITS 20 15 0
__author__ = "Nick Pingel"
__version__ = "1.0"
__email__ = "Nickolas.Pingel@anu.edu.au
__status__ = "Production"
"""

## imports
from astropy.io import fits 
from astropy import units as u 
from astropy.convolution import Gaussian1DKernel
from astropy.convolution import Box1DKernel
from spectral_cube import SpectralCube
import numpy as np
from SPSModule import computePS
## set numpy errors to ignore
np.seterr(divide='ignore', invalid='ignore')
import sys
import pickle
import os
import matplotlib.pyplot as pyplot
from matplotlib.patches import Ellipse
import matplotlib 
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

matplotlib.rc('font', family='sans-serif') 
matplotlib.rc('font', serif='Helvetica Neue') 
matplotlib.rc('text', usetex='false') 
matplotlib.rc('xtick.major.width')
#matplotlib.rcParams['contour.negative_linestyle']= 'solid'
#matplotlib.rcParams.update({'font.size': 14})
# These are the "Tableau 20" colors as RGB.
tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]

# Scale the RGB values to the [0, 1] range, which is the format matplotlib 
# accepts.
for i in range(len(tableau20)):
    r, g, b = tableau20[i]
    tableau20[i] = (r / 255., g / 255., b / 255.)

## function to propagate the individual SPS uncertainties for uncertainty on weighted average value
def propagateError(slopes, errors):
	## compute sum of errors
	err_factor = 1.0/np.sum(errors)
	## compute the numerator
	num_val = 0
	## if computing weighted sum, propagate errors
	if len(slopes) > 1:
		#for n in range(0, len(slopes)):
		#	slope_val = slopes[n]
		#	err_val = errors[n]
		#	num_val += (slope_val**2 * err_val**2)
		## compute final uncertainty value 
		#finalErr_val = err_factor * np.sqrt(num_val)
		finalErr_val = 1/np.sqrt(np.sum(1/errors**2))
		return finalErr_val
	else:
		## if only one value, then no propagation to be done
		return errors[0]

### function to down-sample spectral axis 
def spectralDownSample(current_res, target_res, fwhm_factor, cube):
	current_res = current_res * u.km/u.s
	target_res = target_res * u.km/u.s
	orig_pixel_scale = cube.header['CDELT3']/1000.*u.km/u.s
	oldVelAxis = cube.spectral_axis/1000.*u.km/u.m
	gaussian_width = ((target_res**2 - current_res**2)**0.5 / orig_pixel_scale / fwhm_factor)
	#width = target_res/ orig_pixel_scale
	newVelAxis = np.linspace(oldVelAxis[0], oldVelAxis[-1], np.int(np.ceil(np.abs(oldVelAxis[-1] - oldVelAxis[0])/target_res)))
	kernel = Gaussian1DKernel(gaussian_width)
	newCube = cube.spectral_smooth(kernel)
	interpCube = newCube.spectral_interpolate(newVelAxis, suppress_smooth_warning = True) ## still missing structure in final integrated image -- use a more precise interpolation method???
	return interpCube

## function for progress bar that informats user              
def progressBar(value, endvalue, vel_res, bar_length=20):
	percent = float(value) / endvalue
	arrow = '-' * int(round(percent * bar_length)-1) + '>'
	spaces = ' ' * (bar_length - len(arrow))
	sys.stdout.write("Percent of slopes calculated at resolution %.2f [km/s]:" % (vel_res) + "[{0}] {1}%\r".format(arrow + spaces, int(round(percent * 100))))
	sys.stdout.flush() 

## function to produce plot
def makePlot(slopeArr, errorArr, xAxis, saveFlag, finalPlotFlag):
	## set plotting parameters
	majorYLocFactor = 10**np.floor(np.log10(np.abs(np.max(slopeArr))))/4.
	minorYLocFactor = majorYLocFactor/10.
	majorXLocFactor = 10**np.floor(np.log10(np.abs(np.max(xAxis))))/4.
	minorXLocFactor = majorXLocFactor/10.

	majorYLocator = MultipleLocator(majorYLocFactor)
	majorYFormatter = FormatStrFormatter('%.1f')
	minorYLocator = MultipleLocator(minorYLocFactor)
	majorXLocator = MultipleLocator(majorXLocFactor)
	majorXFormatter = FormatStrFormatter('%.1f')
	minorXLocator = MultipleLocator(minorXLocFactor)
	## do the plotting
	fig, ax = pyplot.subplots()
	pyplot.errorbar(xAxis, (-1)*np.array(slopeArr), yerr=errorArr, fmt='o', linewidth = 2, color = tableau20[0])

	ax.yaxis.set_major_locator(majorYLocator)
	ax.yaxis.set_major_formatter(majorYFormatter)
	ax.yaxis.set_minor_locator(minorYLocator)
	ax.xaxis.set_major_locator(majorXLocator)
	ax.xaxis.set_major_formatter(majorXFormatter)
	ax.xaxis.set_minor_locator(minorXLocator)
	pyplot.ylabel(r'<$\gamma$>')
	## if output is desired, save plot
	if saveFlag == 1 and finalPlotFlag == 1:
		pyplot.xlabel(r'$\Delta v$ [km/s]')
		pyplot.savefig('VCA.pdf')
	elif saveFlag == 1 and finalPlotFlag == 0:
		pyplot.xlabel(r'Velocity [km/s]')
		pyplot.savefig('SPS_vs_AllChans.pdf')
	elif saveFlag == 0 and finalPlotFlag == 1:
		pyplot.xlabel(r'$\Delta v$ [km/s]')
	else:
		pyplot.xlabel(r'Velocity [km/s]')
	pyplot.show()
	pyplot.clf()
	pyplot.close()



## get command line arguments
fitsName = sys.argv[1] ## path to FITS file

## total number of VCA bins
N_VCA = np.int(sys.argv[2])

## total number of spatial freq bins (or elliptical annuli)
N_bins = np.int(sys.argv[3]) 

## save plot flag
saveFlag = np.int(sys.argv[4])

## read in the FITS file and collect information from header like the size, pixel resolution, and beam resolution
hdu = fits.open(fitsName)

## get spatial dimension size
raSize = hdu[0].header['NAXIS1']
decSize = hdu[0].header['NAXIS2']
## get spectral dimension size 
specSize = hdu[0].header['NAXIS3']
## get spectral axis resolution
velRes = hdu[0].header['CDELT3']/1000. ## units of km/s
## get reference velocity & pixel
refVel = hdu[0].header['CRVAL3']/1000. ## units of km/s
refVelPix = hdu[0].header['CRPIX3']
## compute initial velocity
initVel = refVel - refVelPix*velRes
lastVel = initVel + velRes * specSize

## get pixel size 
pixRes = np.abs(hdu[0].header['CDELT1'])

# get beam resolution (major axis)
angRes = hdu[0].header['BMAJ']

## compute the velocity resolutions based on user provided bins
totVel = np.abs(lastVel - initVel)
deltaV = np.abs(lastVel - initVel) / (N_VCA - 1)
print('The resolution of the spectral axis will be down-sampled from %.2f [km/s] to %.2f [km/s] in %s intervals of %.2f [km/s]' % (np.abs(velRes), totVel, np.str(N_VCA), deltaV))
channelWidthArr = np.linspace(np.abs(velRes), totVel, N_VCA) 

## re-read input fits as a spectral-cube object
hdu.close()
origCube = SpectralCube.read(fitsName)

## specify global Gaussian factors
fwhm_factor = np.sqrt(8 * np.log(2))

## create lists to store the weighted average slope values
finalSlopes = []
finalErrors = []

cnt = 0 
## Now, loop over each velocity resolution to down-sample the velocity axis and compute the SPS at each velocity resolution 
for i in range(0, N_VCA):
	## if first iteration, we're measuring the SPS slope at each velocity channel at native resolution
	## As such, make a plot to see how the SPS slope varies with velocity channel before returning the weighted (by uncertainty) average value 
	if i == 0:
		## extract data at native velocity resolution
		origData = origCube.unmasked_data[:,:,:,]

		## create list to store SPS slope & uncertainty values
		slopeList = []
		errorList = []

		## loop through each velocity channel and compute the SPS
		for z in range(0, origData.shape[0]):
			## find and replace NaN values
			intImage = np.copy(origData[z])
			where_are_NaNs = np.isnan(intImage)
			intImage[where_are_NaNs] = 0.
			slope, error = computePS(intImage, N_bins, pixRes, angRes)
			slopeList.append(slope)
			errorList.append(error)

			## provide update via a progress bar:
			if z % 10 == 0:
				progressBar(z, origData.shape[0], channelWidthArr[i])
		## constuct velocity axis
		velocityAxis = np.linspace(initVel, lastVel, specSize)

		## make plot here
		makePlot(slopeList, errorList, velocityAxis, saveFlag, 0)

		## if output desired (savePlot == 1), save SPS slopes vs. native velocity channels
		if saveFlag == 1:
			saveFile = 'VCA.pkl'
			with open(saveFile, "wb") as f:
				pickle.dump([slopeList, errorList, velocityAxis], f)


		## cast into arrays
		slopeArr = np.array(slopeList)
		errorArr = np.array(errorList)

		## compute weighted average
		aveSlope = np.average(slopeArr, weights = 1/errorArr**2)
		## propagate the uncertainty
		propErr = propagateError(slopeArr, errorArr)

		## append values to lists
		finalSlopes.append(aveSlope)
		finalErrors.append(propErr)

		## reset for next iteration
		del slopeList
		del errorList
		del slopeArr
		del errorArr

	## if we're not at the first iteration, then first down-sample the spectral axis 
	else:
		current_res = channelWidthArr[0]
		target_res = channelWidthArr[i]
		## down-sample here
		print('Down-sampling spectral axis from %.2f [km/s] to %.2f [km/s]' %  (current_res, target_res))
		newCube = spectralDownSample(current_res, target_res, fwhm_factor, origCube)

		## extract data at new velocity resolution
		newData = newCube.unmasked_data[:, :, :]
		#if i % 5 == 0:
			#fits.writeto('test_%s.fits' % np.str(i), np.array(newData))

		## create list to store SPS slope & uncertainty values
		slopeList = []
		errorList = []
		## loop through each velocity channel and compute the SPS
		for z in range(0, newData.shape[0]):
			## find and replace NaN values
			intImage = np.copy(newData[z])
			where_are_NaNs = np.isnan(intImage)
			intImage[where_are_NaNs] = 0.
			slope, error = computePS(intImage, N_bins, pixRes, angRes)
			slopeList.append(slope)
			errorList.append(error)
			## provide update via a progress bar:
			if z % 10 == 0:
				progressBar(z, newData.shape[0], target_res)

		## cast into arrays
		slopeArr = np.array(slopeList)
		errorArr = np.array(errorList)

		## compute weighted average
		aveSlope = np.average(slopeArr, weights = 1/errorArr**2)
		## propagate the uncertainty
		propErr = propagateError(slopeArr, errorArr)

		## append values to lists
		finalSlopes.append(aveSlope)
		finalErrors.append(propErr)

		## reset for next iteration
		del slopeList
		del errorList
		del slopeArr
		del errorArr
		cnt+=1

## make final plot here
print('trying to plot')
makePlot(finalSlopes, finalErrors, channelWidthArr, saveFlag, 1)
print('made plot')

## if output is desired (saveFlag == 1), save out final VCA variables 
if saveFlag == 1:
	print('trying to save flag')
	saveFile = 'VCA.pkl'
	with open(saveFile, "wb") as f:
		pickle.dump([finalSlopes, finalErrors, channelWidthArr], f)
	print('trying to save flag')








