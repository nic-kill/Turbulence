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

name='smc_grid_x2_y2'
#if not simulating specify the location of the .fits to source for the vca
sourcefits = '/avatar/nickill/smc/grid_cubes/'+name+'.fits'
figtitle = name
figsaveloc = '/priv/myrtle1/gaskap/nickill/smc/vca/turbustatoutput/' + name + 'witherrors'
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
#DO VCA
#set arrays for VCA calc
vca_array=np.zeros(3)
#arlen=int(len(cube.data[:,0,0]))/10

#do full thickness mom0 SPS and add to array first
#import data and compute moment 0
moment0=cube.moment(order=0)

#compute SPS, add in distance at some point as parameter
pspec = PowerSpectrum(moment0)
pspec.run(verbose=False, xunit=u.pix**-1)
vca_array=np.vstack((vca_array,[pspec.slope,len(cube[:,0,0]),pspec.slope_err]))

#iterate VCA over fractions of the total width of the PPV cube
for i in [128,64,32,16,8,4,2,1]:
    cube.allow_huge_operations=True
    downsamp_cube = cube.downsample_axis(i, axis=0)
    downsamp_cube.allow_huge_operations=True
    vca = VCA(downsamp_cube)
    vca.run(verbose=False, beam_correct=correctbeam)
    vca_array=np.vstack((vca_array,[vca.slope,i,vca.slope_err]))
vca_array=vca_array[1:,:]

#save the array for future plotting without recomputing
np.save(figsaveloc, vca_array)


#for i in [128,64,32,16,8,4,2,1]:
#    vca = VCA(cube, channel_width=i)
#    vca.run(verbose=False, beam_correct=correctbeam)
#    vca_array=np.vstack((vca_array,[vca.slope,i]))
#vca_array=vca_array[1:,:]

#for i in [0.1*arlen,0.2*arlen,0.3*arlen,0.4*arlen,0.5*arlen,0.6*arlen,0.7*arlen,0.8*arlen,0.9*arlen,1*arlen]:
#    vca = VCA(cube, channel_width=i * u.km / u.s)
#    vca.run(verbose=False, beam_correct=correctbeam)
#    vca_array=np.vstack((vca_array,[vca.slope,i]))
#vca_array=vca_array[1:,:]
##########################

##########################
#Plot the VCA by channel

specin=vca_array[:,0]
thickness=vca_array[:,1]
upper_err=vca_array[:,0]+vca_array[:,2]
lower_err=vca_array[:,0]-vca_array[:,2]

plt.title('VCA:' + name)
plt.xlabel('Channel Thickness (km/s)')
plt.ylabel('Spectral Index')
plt.tight_layout()
plt.fill_between(thickness, upper_err, lower_err)
plt.scatter(thickness,specin)
plt.tight_layout()
plt.savefig(figsaveloc + '.png', format='png')