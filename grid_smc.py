import glob
import os
import subprocess,shlex,shutil
import sys
from astropy.io import fits
from spectral_cube import SpectralCube
import numpy as np

splitfactor=7

wholecube=fits.open('/avatar/naomi/ASKAP/SMC/SB_8906/SMC_8906.lsr.K.fits')
print(wholecube[0].shape)

xlen=len(wholecube[0].data[0,:,0])
ylen=len(wholecube[0].data[0,0,:])

for i in np.arange(splitfactor+1):
	xax=np.append(xax,i*xlen/splitfactor)

for i in np.arange(splitfactor+1):
       	yax=np.append(yax,i*ylen/splitfactor)

#yax=[0,ylen/3,(ylen/3)*2,ylen]

wholecube.close()
wholespeccube=SpectralCube.read('/avatar/naomi/ASKAP/SMC/SB_8906/SMC_8906.lsr.K.fits')

for j in np.arange(0,splitfactor):
	for i in np.arange(0,splitfactor):
		sub=wholespeccube.subcube(xlo=int(xax[i]), xhi=int(xax[i+1]), ylo=int(yax[j]), yhi=int(yax[j+1]), zlo=55, zhi=310, rest_value=None)
		sub.write('/avatar/nickill/smc/grid_cubes/smc_grid_x'+str(i)+'_y'+str(j)+'.fits')
		print('done x'+str(i)+' y'+str(j))









