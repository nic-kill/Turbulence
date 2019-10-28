#!/bin/csh

#miriad task to crop in space and velocity

cd /priv/myrtle1/gaskap/nickill/smc/
fits in='/priv/myrtle1/gaskap/GASSIII/GASS_FULL_FK5_ZEA.fits' out=miriadcrop.tmp.image op=xyin
imsub region='images(685,1110)' in='miriadcrop.tmp.image' out=miriadcrop2.tmp.image 
moment in=miriadcrop2.tmp.image mom=0 out=miriadcrop3.tmp.image
fits in='miriadcrop3.tmp.image' out='testcrop.fits' op=xyout line=velocity
rm -rf miriadcrop.tmp.image
rm -rf miriadcrop2.tmp.image
rm -rf miriadcrop3.tmp.image











