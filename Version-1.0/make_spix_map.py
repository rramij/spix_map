## Written by Ramij Raja
## Date : 27 March 2024
###########################
################### Main Code #################################################
################################################################################
# Masking image pixels
mimage = []
for i in range(len(image)):
  data, hdr = getdata(image[i], header=True) # read data
  mask_level = sig * rms[i] # calculate mask level below which pixels will be NaN
  data[data < mask_level] = np.nan # make NaN pixels
  # Mask image beyond the central part
  x_lim = [int(hdr['CRPIX1']-(imsize/2)), int(hdr['CRPIX1']+(imsize/2))]
  y_lim = [int(hdr['CRPIX2']-(imsize/2)), int(hdr['CRPIX2']+(imsize/2))]
  a = data[0,0,0:y_lim[0],:]
  a[a != np.nan] = np.nan
  b = data[0,0,y_lim[1]:,:]
  b[b != np.nan] = np.nan
  c = data[0,0,:,0:x_lim[0]]
  c[c != np.nan] = np.nan
  d = data[0,0,:,x_lim[1]:]
  d[d != np.nan] = np.nan
  # Write FITS image
  fits.writeto('image'+str(i)+'.fits', data, hdr, overwrite=True) 
  mimage.append('image'+str(i)+'.fits')


# Spectral Index Map making
spix_im = []
spix_err_im = []

for i in range(len(mimage)-1):
  # read data
  data0, hdr0 = getdata(mimage[i], header=True)
  data1, hdr1 = getdata(mimage[i+1], header=True)
  # calculate spix
  alpha = np.log(data0/data1) / np.log(freq[i]/freq[i+1])
  # mask +ve alpha values
  alpha[alpha >= 0] = np.nan
  fits.writeto(pre+'_spix_map_'+str(int(freq[0]))+'_'+str(int(freq[1]))+'MHz_'+suf+'.fits',
 alpha, hdr0, overwrite=True)
  spix_im.append(pre+'_spix_map_'+str(int(freq[0]))+'_'+str(int(freq[1]))+'MHz_'+suf+'.fits')
  # calculate spix err
  delta_S0 = np.sqrt( rms[i]**2 + (f_err[i] * data0)**2 )
  delta_S1 = np.sqrt( rms[i+1]**2 + (f_err[i+1] * data1)**2 )
  delta_alpha = (1.0/np.log(freq[i]/freq[i+1])) * np.sqrt( (delta_S0/data0)**2 + (delta_S1/data1)**2 )
  # mask +ve alpha value pixels
  mask = np.isnan(alpha)
  delta_alpha[mask] = np.nan
  fits.writeto(pre+'_spix_err_map_'+str(int(freq[0]))+'_'+str(int(freq[1]))+'MHz_'+suf+'.fits',
   delta_alpha, hdr0, overwrite=True)
  spix_err_im.append(pre+'_spix_err_map_'+str(int(freq[0]))+'_'+str(int(freq[1]))+'MHz_'+suf+'.fits')



# Later add curvature CODE here as well



print(' ')
print('###########################   Spectral Index/Curvature maps created   #########################')
# End




