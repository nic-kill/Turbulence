#!/bin/csh

cd /avatar/nickill/smc/diagnostic_cubes/

rm -rf /avatar/nickill/smc/diagnostic_cubes/smc_masked.im
rm -rf /avatar/nickill/smc/diagnostic_cubes/smctmpim

fits in='/avatar/naomi/ASKAP/SMC/SB_8906/SMC_8906.lsr.K.fits' out='/avatar/nickill/smc/diagnostic_cubes/smctmpim' op=xyin
fits in='/priv/myrtle1/gaskap/nickill/smc/diagnostic_cubes/smc_gain_mask_600chan.fits' out='/avatar/nickill/smc/diagnostic_cubes/gaintmpim' op=xyin
maths exp='<smctmpim>*1' mask='(<gaintmpim>.gt.7e-2)' out='/avatar/nickill/smc/diagnostic_cubes/smc_masked.im'
fits in='/avatar/nickill/smc/diagnostic_cubes/smc_masked.im' out='/avatar/nickill/smc/diagnostic_cubes/smc_masked_0.07.fits' op=xyout

rm -rf /avatar/nickill/smc/diagnostic_cubes/gaintmpim
rm -rf /avatar/nickill/smc/diagnostic_cubes/smc_masked.im
rm -rf /avatar/nickill/smc/diagnostic_cubes/smctmpim
















