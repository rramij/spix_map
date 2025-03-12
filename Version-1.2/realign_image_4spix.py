# To re-align images for spix maps and others
# Written by Ramij Raja
# Date: 5 May 2024
# Run on Conda2 or 3; required:numpy version > 1.12.0
#########################################
#############################
# Main Code 
#############################
for i in range(len(image)-1):
  data, hdr = getdata(image[i+1], header=True)
  # Shift image in +x to shift right and +y to shift up
  # axis = 3 is x shift; axis = 2 is y shift
  data1 = np.roll(data, (int(round(xy_corr[i][0])),int(round(xy_corr[i][1]))), axis=(3,2))
  fits.writeto(image[i+1].replace(".fits", "_realign.fits"), data1, hdr, overwrite=True)
  print('Aligned Image: '+image[i+1].replace(".fits", "_realign.fits"))
  image[i+1] = image[i+1].replace(".fits", "_realign.fits")


print(' ')
print('##  Images Aligned  ##')
# THE END






