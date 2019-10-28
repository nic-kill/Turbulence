#!/bin/csh

#miriad task to crop in space and velocity
#format = boxes(xmin,ymin,xmax,ymax)(z1,z2)
#format = images(z1,z2) Select image planes z1 to z2 inclusive. z2 defaults to z1.


cd /priv/myrtle1/gaskap/nickill/smc/
fits in='/avatar/naomi/ASKAP/SMC/SB_8906/SMC_8906.lsr.K.fits' out=miriadcrop.tmp.image op=xyin
imsub region='boxes(1972,2349,2788,3285)(55,309)' in='miriadcrop.tmp.image' out=miriadcrop2.tmp.image 
fits in='miriadcrop2.tmp.image' out='smc_askap_northbar.fits' op=xyout
rm -rf miriadcrop.tmp.image
rm -rf miriadcrop2.tmp.image









