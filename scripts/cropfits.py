#script to run in CASA that will clip in spatial and then frequency space exporting back into .fits in velocity space
#added drodeg=true and history=False in export to try and fix VCA issue
print(selfits)
print(outfits)
ia.open(selfits)
cropped = ia.crop(region='centerbox[[1685pix, 2073pix], [200pix, 200pix]]', axes=[])
cropped2 = cropped.crop(outfile="cropped"+outfits+".image", chans='48~80', axes=[])
exportfits(imagename='cropped'+outfits+'.image', fitsimage=outfits+'cropped.fits',velocity=True,dropdeg=True,history=False)
