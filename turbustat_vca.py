#MAKE PPV SIM CUBE
import matplotlib.pyplot as plt
from astropy.io import fits
import turbustat
from turbustat.simulator import make_3dfield, make_ppv
import astropy.units as u
from spectral_cube import SpectralCube
from turbustat.io.sim_tools import create_image_header, create_cube_header
from radio_beam import Beam
import numpy as np
from turbustat.statistics import VCA, PowerSpectrum
from turbustat.statistics.apodizing_kernels import CosineBellWindow, TukeyWindow, HanningWindow, SplitCosineBellWindow

#########################
####SET PARAMETERS
dosim = False
simsaveloc = '/priv/myrtle1/gaskap/nickill/smc/vca/simcubes/redo_vel2.4_dens3_gaussiankernel4.fits'

addbeamhead = False
addnoise=False
#correct beam below doesn't work for sim data, only has worked on real data where i have added in the  beam with addbeamhead so far
correctbeam=False
#taper is not yet working can't overwite the data when read from spectral cube, may need to use fits.open and check speccube docs on large datasets
taper=False

name='smc_grid_x0_y0'
#if not simulating specify the location of the .fits to source for the vca
sourcefits = '/avatar/nickill/smc/grid_cubes/'+name+'.fits'
figtitle = name
figsaveloc = '/priv/myrtle1/gaskap/nickill/smc/vca/turbustatoutput/' + name
#########################



#########################
######Simulation
if dosim == True:
	velocity = make_3dfield(100, powerlaw=2.4, amp=1., randomseed=98734) * u.km / u.s
	density = make_3dfield(100, powerlaw=3, amp=1., randomseed=328764) * u.cm**-3
	
	#delete negative denisties by shifitng up a standard dev and zeroing <0 values
	density += density.std()
	density[density.value < 0.] = 0. * u.cm**-3

	#below sets the vel range and chan thickness
	cube_hdu = make_ppv(velocity, density, los_axis=0, T=100 * u.K, chan_width=0.1 * u.km / u.s, v_min=-10 * u.km / u.s, v_max=10 * u.km / u.s)

	#read fits file with SpectralCube
	cube = SpectralCube.read(cube_hdu)	

	#adds beam info to header
	if addbeamhead==True:
		#specify a beam size
		madebeam=Beam(30*u.arcsec)
		#add in beam
		cube.with_beam(madebeam)

if dosim==True:
	cube.write(simsaveloc)
#########################


########################
#PRE-PROCESSING

if dosim==False:
	#cube = fits.open(sourcefits)[0]
	cube = SpectralCube.read(sourcefits)

if dosim==True:
        cube = fits.open(simsaveloc)[0]

if addbeamhead==True:
        #specify a beam size
        madebeam=Beam(30*u.arcsec)
        cube=SpectralCube.read(cube)
	#add in beam
        cube.with_beam(madebeam)

if addnoise==True:
	sigma = 0.1
	cube.data = cube.data + np.random.normal(0., sigma, size=cube.shape)
	figsaveloc = '/priv/myrtle1/gaskap/nickill/smc/vca/turbustatoutput/' + name + '_noise.png'

if taper==True:
	#unmasked data has shape (spectral,n_y,n_x)
	shape = (len(cube.unmasked_data[0,:,0]), len(cube.unmasked_data[0,0,:]))
	taper = TukeyWindow(alpha=0.6)
	window = taper(shape)
	for slice in range(len(cube.unmasked_data[:,0,0])):
		cube.unmasked_data[slice,:,:] = cube.unmasked_data[slice,:,:]*(1-window)
########################



########################
def do_vca(vcacube, array_save_loc, fig_save_loc):
    """This function greets to
    the person passed in as
    parameter"""
    
    vca_array=np.zeros(3)
    #arlen=int(len(vcacube.data[:,0,0]))/10
    
    #do full thickness mom0 SPS and add to array first
    #import data and compute moment 0
    moment0=vcacube.moment(order=0)
    
    #compute SPS, add in distance at some point as parameter
    pspec = PowerSpectrum(moment0)
    pspec.run(verbose=False, xunit=u.pix**-1)
    vca_array=np.vstack((vca_array,[pspec.slope,len(vcacube[:,0,0]),pspec.slope_err]))

    #iterate VCA over fractions of the total width of the PPV vcacube
    for i in [128,64,32,16,8,4,2,1]:
        vcacube.allow_huge_operations=True
        downsamp_vcacube = vcacube.downsample_axis(i, axis=0)
        downsamp_vcacube.allow_huge_operations=True
        vca = VCA(downsamp_vcacube)
        vca.run(verbose=False, beam_correct=correctbeam, save_name=fig_save_loc+'_thickness'+str(i)+'.png')
        vca_array=np.vstack((vca_array,[vca.slope,i,vca.slope_err]))
    vca_array=vca_array[1:,:]

    #save the array for future plotting without recomputing
    np.save(array_save_loc, vca_array)

for j in np.arange(0,7):
	for i in np.arange(0,7):
		cube = SpectralCube.read('/avatar/nickill/smc/grid_cubes/smc_grid7x7_masked_x'+str(i)+'_y'+str(j)+'.fits')
		arrayloc = '/priv/myrtle1/gaskap/nickill/smc/vca/turbustatoutput/smc_grid7x7_masked_x'+str(i)+'_y'+str(j)
		figloc = '/priv/myrtle1/gaskap/nickill/smc/vca/turbustatoutput/smc_grid7x7_masked_x'+str(i)+'_y'+str(j)
		do_vca(cube,arrayloc,figloc)                
		print('done x'+str(i)+' y'+str(j))




