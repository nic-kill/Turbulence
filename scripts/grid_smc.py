import glob
import os
import subprocess,shlex,shutil
import sys
from astropy.io import fits
from spectral_cube import SpectralCube
import numpy as np

################
#Parameters
#define the number of subcubes per axis
splitfactor=31

#specify source cube location
sourcefile='/avatar/nickill/smc/diagnostic_cubes/smc_masked_0.07.fits'

cubenameprefix='/avatar/nickill/smc/grid_cubes/smc_grid31x31_masked'
###################



###################
##Find dimensions

wholecube=fits.open(sourcefile)
print(wholecube[0].shape)

xlen=len(wholecube[0].data[0,:,0])
ylen=len(wholecube[0].data[0,0,:])

xax=[]
for i in np.arange(splitfactor+1):
	xax=np.append(xax,i*xlen/splitfactor)

yax=[]
for i in np.arange(splitfactor+1):
       	yax=np.append(yax,i*ylen/splitfactor)

#yax=[0,ylen/3,(ylen/3)*2,ylen]

wholecube.close()
##################


##################
#Make mom0 to overlay regions on and split off subregions
wholecube=SpectralCube.read(sourcefile)

#make the mom0 to overwrite
moment0=wholecube.moment(order=0)

for j in np.arange(0,splitfactor):
	for i in np.arange(0,splitfactor):
		print('starting x'+str(i)+' y'+str(j))
		#overwrite region boundaries with really high values
		moment0.array[:,int(xax[i])]=99999999
		moment0.array[int(yax[j]),:]=99999999
		#split off sub regions
		sub=wholecube.subcube(xlo=int(xax[i]), xhi=int(xax[i+1]), ylo=int(yax[j]), yhi=int(yax[j+1]), zlo=55, zhi=310, rest_value=None)
		sub.write(cubenameprefix+'_x'+str(i)+'_y'+str(j)+'.fits')
		print('done x'+str(i)+' y'+str(j))
moment0.write(cubenameprefix+'_regionoutlines.fits')
#################







