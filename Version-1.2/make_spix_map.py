## Written by Ramij Raja
## Date : 11 March 2025
###########################
################### Main Code #################################################
################################################################################
# Masking image pixels
print(' ')
print('Masking image pixels beyond central region...')
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
print(' ')
print('Starting spix map making...')
ref_hdr = getheader(mimage[0])

spix_im = []
spix_err_im = []

for i in range(len(mimage)-1):
  # read data
  data0, hdr0 = getdata(mimage[i], header=True)
  data1, hdr1 = getdata(mimage[i+1], header=True)
  # calculate spix
  alpha = np.log(data0/data1) / np.log(freq[i]/freq[i+1])
  # mask +ve/-ve alpha values
  if maskrange == True:
    alpha[alpha >= max(maskpix)] = np.nan
    alpha[alpha <= min(maskpix)] = np.nan
  fits.writeto(pre+'_spix_map_'+str(int(freq[i]))+'_'+str(int(freq[i+1]))+'MHz_'+suf+'.fits',
 alpha, ref_hdr, overwrite=True)
  spix_im.append(pre+'_spix_map_'+str(int(freq[i]))+'_'+str(int(freq[i+1]))+'MHz_'+suf+'.fits')
  print(' ')
  print('Spix map: '+pre+'_spix_map_'+str(int(freq[i]))+'_'+str(int(freq[i+1]))+'MHz_'+suf+'.fits')
  # calculate spix err
  delta_S0 = np.sqrt( rms[i]**2 + (f_err[i] * data0)**2 )
  delta_S1 = np.sqrt( rms[i+1]**2 + (f_err[i+1] * data1)**2 )
  delta_alpha = (1.0/abs(np.log(freq[i]/freq[i+1]))) * np.sqrt( (delta_S0/data0)**2 + (delta_S1/data1)**2 )
  # match spix map mask pixels
  mask = np.isnan(alpha)
  delta_alpha[mask] = np.nan
  fits.writeto(pre+'_spix_err_map_'+str(int(freq[i]))+'_'+str(int(freq[i+1]))+'MHz_'+suf+'.fits',
   delta_alpha, ref_hdr, overwrite=True)
  spix_err_im.append(pre+'_spix_err_map_'+str(int(freq[i]))+'_'+str(int(freq[i+1]))+'MHz_'+suf+'.fits')
  print('Spix Err map: '+pre+'_spix_err_map_'+str(int(freq[i]))+'_'+str(int(freq[i+1]))+'MHz_'+suf+'.fits')

print(' ')
print('Spix maps DONE! Starting Carvature maps...')
# Curvature CODE 
for i in range(len(spix_im)-1):
  # read data
  spix0 = getdata(spix_im[i])
  spix1 = getdata(spix_im[i+1])
  # calculate curvature
  curv = spix1 - spix0
  fits.writeto(pre+'_curve_map_'+str(int(freq[i]))+'_'+str(int(freq[i+2]))+'MHz_'+suf+'.fits', curv, ref_hdr, overwrite=True)
  print(' ')
  print('Curvature Map: '+pre+'_curve_map_'+str(int(freq[i]))+'_'+str(int(freq[i+2]))+'MHz_'+suf+'.fits')
  # calculate curvature
  a = - np.log(freq[i+1]/freq[i+2]) * rms[i]
  b = np.log(freq[i]/freq[i+2]) * rms[i+1]
  c = - np.log(freq[i]/freq[i+1]) * rms[i+2]
  d = np.log(freq[i]/freq[i+1]) * np.log(freq[i+1]/freq[i+2])
  im0 = getdata(mimage[i])
  im1 = getdata(mimage[i+1])
  im2 = getdata(mimage[i+2])
  curv_e = np.sqrt((a/im0)**2 + (b/im1)**2 + (c/im2)**2)/d
  # mask match curv pixels
  mask = np.isnan(curv)
  curv_e[mask] = np.nan
  fits.writeto(pre+'_curve_err_map_'+str(int(freq[i]))+'_'+str(int(freq[i+2]))+'MHz_'+suf+'.fits', curv_e, ref_hdr, overwrite=True)
  print('Curvature Err Map: '+pre+'_curve_err_map_'+str(int(freq[i]))+'_'+str(int(freq[i+2]))+'MHz_'+suf+'.fits')


print(' ')
print('##   Spectral Index/Curvature maps created   ##')
# End




