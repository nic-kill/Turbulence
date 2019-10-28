#!/bin/csh

cd /avatar/nickill/smc/
mossen in='/avatar/naomi/ASKAP/SMC/SB_8906/SMC_8906_00.cln' sen='/avatar/nickill/smc/grid_cubes/smc_sen.image' gain='/avatar/nickill/smc/grid_cubes/smc_gain.image'
fits in='/avatar/nickill/smc/grid_cubes/smc_sen.image' out='/avatar/nickill/smc/grid_cubes/smc_sen.fits' op=xyout
fits in='/avatar/nickill/smc/grid_cubes/smc_gain.image' out='/avatar/nickill/smc/grid_cubes/smc_gain.fits' op=xyout

rm -rf /avatar/nickill/smc/grid_cubes/smc_sen.image
rm -rf /avatar/nickill/smc/grid_cubes/smc_gain.image
