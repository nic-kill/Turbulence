"""
origData[z], N_bins, raSize, decSize, pixRes, angRes
5/31/19
Modular version of SPS_2D.py. That is --- provided a 2D image, the number of annuli, pixel-resolution, and angular resolution --- this script will produce 
a fit to the 2D spatial power spectrum (SPS) without printing any results to terminal nor generating any plots. The SPS is computed by: (1) taking the Fourier 
Transform of the input image and squaring (2) measuring the median power within equal elliptical logarithmic bins of spatial frequency as we're expecting a power law relationatioship (3) and fitting a linear least-squared fit to these power values over 
a logarithmic scale of linear lengths. The slope and uncertainty is returned to the parent program.
Inputs:
(1) - 2D image
(2) - number of bins (annuli)
(3) - pixel resolution [deg]
(4) -  angular resolution [deg]
Usage:
ipython SPS_2D.py /path/to/FITS/image/ 12 1
Example:
ipython SPS_2D.py SMC_ASKAPsoft_JD_2000min_5maj_Parkes_K_MOM0_ATCAREGRID.FITS 15 0
__author__ = "Nick Pingel"
__version__ = "1.0"
__email__ = "Nickolas.Pingel@anu.edu.au
__status__ = "Production"
"""

## imports
from astropy.io import fits 
import scipy.io
import numpy as np
np.seterr(divide='ignore', invalid='ignore')
import sys
import pickle
from scipy.optimize import curve_fit
from scipy import optimize
from scipy.interpolate import *
import matplotlib.pyplot as pyplot
from matplotlib.patches import Ellipse
import matplotlib 
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

matplotlib.rc('font', family='sans-serif') 
matplotlib.rc('font', serif='Helvetica Neue') 
matplotlib.rc('text', usetex='false') 
matplotlib.rc('xtick.major.width')
matplotlib.rcParams['contour.negative_linestyle']= 'solid'
matplotlib.rcParams.update({'font.size': 14})
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

## get command line arguments
fitsName = sys.argv[1] ## path to FITS file

## total number of spatial freq bins (or elliptical annuli)
numAnnuli = np.int(sys.argv[3]) 

"""
method to compute radius of ellipse givien input semi-major/minor axes and angle
"""
def computeEllRadius(maxA, maxB, theta):
	ellRad = maxA*maxB / np.sqrt((maxB*np.cos(theta))**2 + (maxA*np.sin(theta))**2)
	return ellRad 

"""
method to compute radius of given pixel coordinates
"""
def computePixRadius(i, j, cenI, cenJ):
	pixRad = np.sqrt((i - cenI)**2 + (j - cenJ)**2)
	return pixRad

"""
method to compute relevent gridding paramters and coordinates; returns semi-major/minor values in modulus image for lower, middle, 
and upper elliptical annuli (in pixels), the maximum spatial frequency (based on resolution), and pixel increment in Fourier space.
Returns:
	semi-major/minor pixel counts - lower, middle, and upper elliptical annuli boundaries
	semi-major/minor effective spatial frequencies - semi-major set by angular resolution limit, while semi-major set by the ratio of minor/major image pixel size
	modulus pixel increment - set by the inverse of maximum anglular scale of input intensity image
	semi-major/minor spatial frequency size of grid - both set by image pixel size
	angular scale - set by the inverse of middle boundary of spatial frequency pixel vector
Inputs:
	pixRes - the angular resolution of single pixel in original image [deg]
	angRes - the angular resolution of the beam [deg]
	numAnnuli - the number of total elliptical annuli whose widths increase equally in logarithmic space
	maxA - the maximum size of the semi-major axis (1/2 the number of pixels contained in the larger spatial dimension)
	maxB - the maximum size of the semi-minor axis (1/2 of number of pixels contained in shorter spatial dimension)
"""
def computeScales(pixRes, angRes, numAnnuli, maxA, maxB):
	totAngScale = np.deg2rad(maxA*pixRes*2)
	## determine the spatial resolution in the Fourier plane (i.e. in the modulus image)
	minSpatialFreq = 1/totAngScale 
	## compute max and increments in log space of semi-major/minor (A/B) axes
	## In the case where the pixels are independent (same size as beam), set the maximum Effective Spatial Freq to be 
	## the size of the image
	if np.round(pixRes, 3) == np.round(angRes, 3):
		maxEffSpatialFreq_Maj = maxA*minSpatialFreq
		maxEffSpatialFreq_Min = maxB*minSpatialFreq	

	## in the case where pixels are not independent of the beam, set maximum effective Spatial Frequency to be based on 
	## major axis of beam
	else:
		maxEffSpatialFreq_Maj = 1/np.deg2rad(angRes)
		maxEffSpatialFreq_Min = maxEffSpatialFreq_Maj*maxB/maxA
	maxGridSpatialFreq_Maj = maxA*minSpatialFreq ## only up to gridded extent
	maxGridSpatialFreq_Min = maxB*minSpatialFreq

	lowerMaj = []
	lowerMin = []
	midMaj = []
	midMin = []
	upperMaj = []
	upperMin = []

	## determine the increments in logarithmic spatial frequency
	logInc_Major = (np.log(maxEffSpatialFreq_Maj) - np.log(minSpatialFreq))/numAnnuli
	logInc_Minor = (np.log(maxEffSpatialFreq_Min) - np.log(minSpatialFreq))/numAnnuli

	## determine lowest spatial frequency bin
	initLower = np.log(minSpatialFreq)

	#print('Angular_min [deg]: %s' % pixRes)
	#print('Modulus Image pix res: %s' % minSpatialFreq)
	#print('k_max (Major Direction) boudary in Grid [lambda]: %s' % maxGridSpatialFreq_Maj)
	#print('k_max (Minor Direction) boudary in Grid [lambda]: %s' % maxGridSpatialFreq_Min)	
	#print('k_max (Major Direction) boudary in largest elliptical annulus [lambda]: %s' % maxEffSpatialFreq_Maj)
	#print('k_max (Minor Direction) boudary in largest elliptical annulus [lambda]: %s' % maxEffSpatialFreq_Min)
	#print('Angular_max [deg]: %s' % np.rad2deg(totAngScale))
	#print('(Major) Log increment: %s' % logInc_Major)
	#print('(Minor) Log increment: %s' % logInc_Minor)
	for i in range(0, int(numAnnuli)):
	    lowerMaj.append(initLower+(i*logInc_Major))
	    lowerMin.append(initLower+(i*logInc_Minor))
	    midMaj.append(lowerMaj[i]+logInc_Major/2.)
	    midMin.append(lowerMin[i]+logInc_Minor/2.)
	    upperMaj.append(lowerMaj[i]+logInc_Major)
	    upperMin.append(lowerMin[i]+logInc_Minor)

	## convert back to pixel space and return; corresponding angular scale will be computed later on...
	lowerMajPixArr = np.exp(np.array(lowerMaj, dtype='float32'))/minSpatialFreq
	lowerMinPixArr = np.exp(np.array(lowerMin, dtype='float32'))/minSpatialFreq
	midMajPixArr = np.exp(np.array(midMaj, dtype='float32'))/minSpatialFreq
	midMinPixArr = np.exp(np.array(midMin, dtype='float32'))/minSpatialFreq
	upperMajPixArr =np.exp(np.array(upperMaj, dtype='float32'))/minSpatialFreq
	upperMinPixArr =np.exp(np.array(upperMin, dtype='float32'))/minSpatialFreq
	angularScale = np.rad2deg(1/(midMajPixArr*minSpatialFreq))
	return lowerMajPixArr, lowerMinPixArr, midMajPixArr, midMinPixArr, upperMajPixArr, upperMinPixArr, maxEffSpatialFreq_Maj, maxEffSpatialFreq_Min, minSpatialFreq, maxGridSpatialFreq_Maj, maxGridSpatialFreq_Min, angularScale

## linear fitting function
def linFunc(x, slope, b):
	return x*slope+b

"""
Method to perform residual bootstraping on linear fit given data points, their associated uncertainties, and original fit to data
Slope will be computed 1000x
"""
def residBootstrap(scaleArr, medianArr, errArr, fitArr):
	## set total number of iterations and initialize arrays to hold slope, error values
	niters = 100
	finSlopeArr = np.zeros([niters])
	finInterArr = np.zeros([niters])

	## compute inital residuals
	initResids = medianArr - fitArr

	## To bootstrap, first begin loop over total iterations
	for r in range(0, niters):
		## now, for each median power value, randomly select an initial residual
		## by generating a random index array with the same size as the residual array
		indArr = np.random.randint(len(initResids), size=len(initResids))

		## add these selected residuals to original fitted power values
		newFitVals = fitArr + initResids[indArr]

		## recompute the fit
		bootCoeffs, bootMatcov = curve_fit(linFunc, scaleArr, newFitVals,[1,1], sigma = errArr, absolute_sigma = True)
		finSlopeArr[r] = bootCoeffs[0]
		finInterArr[r] = bootCoeffs[1]

	## take final slope value to be the mean of the bootstrapped slopes
	finSlope = np.mean(finSlopeArr)
	finInter = np.mean(finInterArr)

	## compute the variance of the returned slope values
	slopeVar = 1.0/niters * np.sum(finSlopeArr**2) - (1/niters * np.sum(finSlopeArr))**2
	return finSlope, finInter, np.sqrt(slopeVar)



"""
Method that computes the spatial power spectrum of a given 2D input image, number of bins (annuli) pixel size, and angular resolution (in degs)
"""
def computePS(intImage, numAnnuli, pixRes, angRes):

	## determine the maximum semi major and minor values based on size of spatial axes
	## recall that numpy arrays are row major -> first index of .shape is # of rows (or y-axis in image), 
	## while second index is # of columns (or x-axis in image)
	if intImage.shape[0] > intImage.shape[1]:
		maxA = intImage.shape[0]/2
		maxB = intImage.shape[1]/2
	else:
		maxA = intImage.shape[1]/2
		maxB = intImage.shape[0]/2	

	## perform the FFT and construct the power spectrum
	modulusImage = np.abs(np.fft.fftshift(np.fft.fft2(intImage, intImage.shape)))**2
	## get pixel boundaries for elliptical annuli
	lowerMajPixArr, lowerMinPixArr, midMajPixArr, midMinPixArr, upperMajPixArr, upperMinPixArr, maxEffSpatialFreq_Maj, maxEffSpatialFreq_Min, pixInc, maxGridSpatialFreq_Maj, maxGridSpatialFreq_Min, angularScale = computeScales(pixRes, angRes, numAnnuli, maxA, maxB)       

	## central pixel is determined from where the maximum value occrus (index of the DC component)
	cenPix = np.where(modulusImage == np.nanmax(modulusImage))
	cenPix_i = cenPix[1]
	cenPix_j = cenPix[0]
	cenPix_i = cenPix_i[0]
	cenPix_j =cenPix_j[0]
	#print('The center pixels are: %s, %s' % (cenPix_i, cenPix_j))




	##lists to hold results
	medianList = []
	meanList = []
	errList = []
	totPointsList = []
	meanPercentDiffList = []
	fittedMeanPercentDiffList = []

	## loop over each annuli to construct the boundaries and compute the statistics compute statistics
	for ring in range(0, numAnnuli):
		lowerA = lowerMajPixArr[ring]
		lowerB = lowerMinPixArr[ring]
		upperA = upperMajPixArr[ring]
		upperB = upperMajPixArr[ring]
		#print('Computing median value within annulus: %s' % (ring+1))    
		colDenDist = []
		for i in range(0, modulusImage.shape[1]):
			for j in range(0, modulusImage.shape[0]):
				## compute angle
				theta = np.arctan2(j - cenPix_j, i - cenPix_i)
				## compute radius of lower boundary
				rLower = computeEllRadius(lowerA, lowerB, theta)
				## compute radius of upper boundary
				rUpper = computeEllRadius(upperA, upperB, theta)
	        	## compute radius of pixel
				radius = computePixRadius(i, j, cenPix_i, cenPix_j)
				if radius < 1 and ring == 0:
					radius = lowerMinPixArr[ring]
					colDenDist.append(modulusImage[cenPix_j,cenPix_i])
					modulusImage[j,i] = np.float('nan')
				if rLower <= radius < rUpper:
					colDenDist.append(modulusImage[j,i])
					modulusImage[j,i] = np.float('nan')

	    ## append relevant quantities
		colDenDistArr = np.array(colDenDist, dtype='float32')
	    ## remove nan's
		colDenDistArr = colDenDistArr[~np.isnan(colDenDistArr)]
		## compute statistics
		meanVal = np.nanmean(colDenDistArr)
		medianVal = np.nanmedian(colDenDistArr)
		medianList.append(medianVal) 
		meanList.append(meanVal)
		totDataPoints = len(colDenDistArr)
		totPointsList.append(totDataPoints)
		## compute uncertainty via MAD
		errList.append(np.nanmedian(np.abs(colDenDistArr-medianList[ring])))
		#print('Total data points: %s' % len(colDenDistArr))  
		#print(colDenDistArr)
		#print('Mean: %s' % meanVal)
		#print('Median: %s' % medianVal)
		#fits.writeto('/Users/npingel/Desktop/test_%s.fits' % np.str(ring), np.log10(modulusImage))

	#print('Fitting spectra...')

	## error propagation
	errArr = np.array(errList, dtype='float32')
	medianArr = np.array(medianList, dtype='float32')
	meanArr = np.array(meanList, dtype='float32')
	errArr_Log = errArr/(medianArr*np.log(10))
	angScaleLog = np.log10(angularScale)
	medianArrLog = np.log10(medianArr)
	## linear fit
	coeffs, matcov = curve_fit(linFunc, angScaleLog, medianArrLog,[1,1], sigma = errArr_Log, absolute_sigma = True)
	#coeffs,matcov = curve_fit(linFunc, angScaleLog,np.log10(medianArr),[1,1])

	## compute final slope & uncertainty by residual bootstrap
	fitArr = coeffs[0]* angScaleLog + coeffs[1]
	slope, inter, slopeErr = residBootstrap(angScaleLog, medianArrLog, errArr_Log, fitArr)

	## compute error
	#error = np.sqrt(np.diag(matcov))
	#slopeErr = error[0]
	
	## DEBUG PLOTS
	"""
	newFitArr = slope* angScaleLog + inter
	fig, ax = pyplot.subplots()
	pyplot.errorbar(angularScale, np.log10(medianArr), yerr=errArr_Log, fmt='o', linewidth = 2, color = tableau20[0])
	pyplot.plot(angularScale, newFitArr, color='black', linewidth=2, linestyle = '--', label=r'Slope: %.2f' % slope +'+/$-$'+'%.2f' % slopeErr)
	ax.set_xlim(np.max(angularScale)+1/2*np.max(angularScale), np.min(angularScale)-1/2*np.min(angularScale))
	ax.set_xscale('log')
	pyplot.ylabel(r'Log$_{10}$(Power)')
	pyplot.xlabel(r'Angular Scale [deg]')
	pyplot.legend(loc=0, fontsize=14)
	pyplot.show()
	pyplot.clf()
	pyplot.close()
	"""
	## END DEBUG


	

	## return fitted slope and associated uncertainty 
	return slope, slopeErr



