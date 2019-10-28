#!/bin/csh

#miriad task to crop in space and velocity

cd /priv/myrtle1/gaskap/nickill/smc/
fits in='/priv/myrtle1/gaskap/nickill/smc/vca/simcubes/simcube_smc_30arcsecbeamheader_30arcsecmiriadconvol.fits' out=miriadcrop.tmp.image op=xyin
convol map=miriadcrop.tmp.image fwhm=30 out=miriadcrop2.tmp.image options=divide sigma=1
fits in='miriadcrop2.tmp.image' out='/priv/myrtle1/gaskap/nickill/smc/vca/simcubes/simcube_smc_30arcsecbeamheader_30arcsecmiriadconvol_deconv.fits' op=xyout line=velocity
rm -rf miriadcrop.tmp.image
rm -rf miriadcrop2.tmp.image












