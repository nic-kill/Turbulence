#!/bin/csh

#miriad task to crop in space and velocity

cd /priv/myrtle1/gaskap/nickill/smc/
fits in='/avatar/nipingel/ASKAP/SMC/data/smc2019/CUBES/CO-Fields/SMC_SB8906_MSCLEAN_Briggs_r1.1_CONTSUB_chanChunk3_beams25-26.FITS' out=miriadcrop.tmp.image op=xyin
imsub region='boxes(1585,1973,1785,2173)(48,80)' in='miriadcrop.tmp.image' out=miriadcrop2.tmp.image 
fits in='miriadcrop2.tmp.image' out='smcCOmiriadcrop.fits' op=xyout line=velocity
rm -rf miriadcrop.tmp.image
rm -rf miriadcrop2.tmp.image












