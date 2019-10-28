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
splitfactor=7

#specify source cube location
sourcefile='/avatar/nickill/smc/diagnostic_cubes/smc_masked_0.07.fits'
#Naomis original smc cube: '/avatar/naomi/ASKAP/SMC/SB_8906/SMC_8906.lsr.K.fits'

cubenameprefix='/avatar/nickill/smc/grid_cubes/smc_grid7x7_masked'


wholecube=fits.open(sourcefile)
print(wholecube[0].shape)
###################



###################
##Find dimensions
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



#################
#make mom0 and overwrite the mom0 with region


##################



##################
#Split off subregions
wholecube=SpectralCube.read(sourcefile)

for j in np.arange(0,splitfactor):
	for i in np.arange(0,splitfactor):
		sub=wholecube.subcube(xlo=int(xax[i]), xhi=int(xax[i+1]), ylo=int(yax[j]), yhi=int(yax[j+1]), zlo=55, zhi=310, rest_value=None)
		sub.write(cubenameprefix+'_x'+str(i)+'_y'+str(j)+'.fits')
		print('done x'+str(i)+' y'+str(j))
#################








