## Written by Ramij Raja
## Date : 27 March 2024
###########################
## Run on Conda2 or 3; required:numpy version > 1.12.0

import os
import shutil
import numpy as np
from astropy.io import fits
from astropy.io.fits import getdata, getheader
#from astropy.io.fits import getheader
############ INPUTS ##############################
# Image FITS filenames
image = ['A85_400MHz_spix_100luvcut_rob-0.5-MFS-image_bm715631_pbcor.fits',
 'A85_700MHz_spix_100luvcut_rob-0.5-MFS-image_bm715631_pbcor.fits','A85_1280MHz_spix_100luvcut_rob-0.5_nt2-MFS-image_bm715631_PBCOR.fits']

file1 = open('align.dat', 'r')

# Re-align image; [x, y] pixels, +x to shift right and +y to shift up
# For multiple images [[1,2],[3,4]] input format
#xy_corr = [[3,-1]]

# Frequencies in MHz
freq = [400., 700., 1280.]

# Image RMS in Jy
rms = [5.0e-5, 2e-5, 7e-6]

# Mask sigma level
sig = 4

# Reference Frequency in MHz for final image header
# Keep it as the first list element ALWAYS
ref = 400.

# Flux density error
f_err = [0.07, 0.06, 0.05]

# Output file prefix and suffix
pre = "A85" # Cluster Name
suf = "B7x5" # Beam resolution

# Output central spix image size; other pixels will be NaN
imsize = 2000 
################### Main Code #################################################
################################################################################

# First Regrid image
#execfile('regrid_image.py')

def realign(image,eim):
  data, hdr = getdata(image, header=True)
  # Shift image in +x to shift right and +y to shift up
  # axis = 3 is x shift; axis = 2 is y shift
  data1 = np.roll(data, (int(round(eim[0])),int(round(eim[1]))), axis=(3,2))
  fits.writeto(image.replace(".fits", "_realign.fits"), data1, hdr, overwrite=True)


#realign(image,eim4)
dat1 = file1.readlines()

xy_corr = []
for i in range(2,len(dat1)-1):
  a = np.array(dat1[1].split(), dtype=float)
  b = np.array(dat1[i].split(), dtype=float)
  eim = a-b
  eim = np.delete(eim, 0)
  xy_corr.append(eim)

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
print('##   THE END   ##')
print('Enjoy!, ... or Not???')
# End



