## Written by Ramij Raja
## Date : 27 March 2024
###########################
## Run on Conda2 or 3; required:numpy version > 1.12.0

import os
import shutil
import numpy as np
from astropy.io import fits
from astropy.io.fits import getdata

############ INPUTS ##############################
# Image FITS filenames
image = ['ddt_400MHz_rob0_uvsub_nt4-MFS-image_bm8p5_pbcor.fits',
 '36_006_700MHz_rob0_uvsub_nt4-MFS-image_bm8p5_pbcor.fits']

# Re-align image; [x, y] pixels, +x to shift right and +y to shift up
# For multiple images [[1,2],[3,4]] input format
xy_corr = [[3,-1]]

# Frequencies in MHz
freq = [400., 675.]

# Image RMS in Jy
rms = [2.0e-5, 1.5e-5]

# Mask sigma level
sig = 4

# Regrid image reference Frequency in MHz 
#ref = 400.

# Flux density error
f_err = [0.07, 0.06]

# Output file prefix and suffix
pre = "Sausage" # Cluster Name
suf = "B8p5" # Beam resolution

# Output central spix image size; other pixels will be NaN
imsize = 2000 
################### Main Code #################################################
################################################################################

# First Regrid image
#execfile('regrid_image.py')

# First re-align the images
exec(open('realign_image_4spix.py').read())

# Now make Spix maps
exec(open('make_spix_map.py').read())


# Remove intermediate files
for i in mimage:
  if os.path.exists(i):
    # Delete the file
    os.remove(i)
  else: continue  

print(' ')
print('###########################   THE END   #########################')
print('Enjoy!, ... or Not???')
# End



