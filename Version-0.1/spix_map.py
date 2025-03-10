## Written by Ramij Raja
## Date : 01 May 2023
###########################
## USE CASA 4.6 Version
import shutil
############ INPUTS ##############################
# Frequencies in MHz

freq = [323., 700., 1280.]

# Image RMS in Jy

#rms = [100e-6, 20e-6, 10e-6]
rms = [150e-6, 20e-6, 10e-6]

# Mask image Frequency in MHz 
ref = 700.

# Mask Level in Jy
level = 348e-6

# Image FITS filenames

image = ['A85_323MHz_spix_100uv29kl_rob0-image_restore_bm9.fits', 'A85_700MHz_spix_100uv29kl_rob0-MFS-image_restore_bm9.fits', 'A85_1.28GHz_spix_100uv29kl_rob0-MFS-image_restore_bm9.fits']

# Output file prefix and suffix

pre = "A85" # Cluster Name
suf = "B9" # Beam resolution

################### Main Code #################################################
################################################################################

# Import FITS files
for i in image:
	importfits(fitsimage = i,  imagename = i.replace(".fits", ""))
	j = image.index(i)
	image[j] = i.replace(".fits", "")

	
# Regrid images
for i in range(len(freq)):
	if ref == freq[i]: ref_image = image[i]
	else: continue
	
for i in image:
	if i == ref_image : continue
	else: imregrid(imagename = i, template = ref_image, output = i + ".regrid")
	j = image.index(i)
	image[j] = i + ".regrid"	

	
# Creat Mask
if '-' in ref_image :
	shutil.move(ref_image, ref_image.replace("-", "QQQ"))
	ref_image = ref_image.replace("-", "QQQ")	
	ia.open(ref_image)
	ia.calcmask(mask = ref_image + ">" + str(level), name = 'refmask')
	ia.done()
	shutil.move(ref_image, ref_image.replace("QQQ", "-"))
	ref_image = ref_image.replace("QQQ", "-")
else :
	ia.open(ref_image)
	ia.calcmask(mask = ref_image + ">" + str(level), name = 'refmask')
	ia.done()


for i in image:
	if i == ref_image : continue	
	else:
		ia.open(i)
		ia.maskhandler('copy',[ref_image + ":refmask", 'mymask' + str(image.index(i))])
		ia.maskhandler('set', 'mymask' + str(image.index(i)))
		ia.done()

# Spectral Index Map
spix_im = []
spix_err_im = []

for i in range(len(image)-1):
	immath(imagename = [image[i], image[i+1]], mode='spix', outfile = pre + "_spix_map_" +  str(int(freq[i])) + "_" + str(int(freq[i+1])) + "MHz_" + suf)
	
	exportfits(imagename = pre + "_spix_map_" +  str(int(freq[i])) + "_" + str(int(freq[i+1])) + "MHz_" + suf, fitsimage = pre + "_spix_map_" +  str(int(freq[i])) + "_" + str(int(freq[i+1])) + "MHz_" + suf + ".fits", overwrite=True)
	
	spix_im.append(pre + "_spix_map_" +  str(int(freq[i])) + "_" + str(int(freq[i+1])) + "MHz_" + suf)
	
#	shutil.rmtree(pre + "_spix_map_" +  str(int(freq[i])) + "_" + str(int(freq[i+1])) + "MHz_" + suf)

	expr ='(1.0/abs(log(' + str(freq[i]) + '/' + str(freq[i+1]) + ')))*sqrt(pow((' + str(rms[i]) + '/IM0),2)+pow((' + str(rms[i+1]) + '/IM1),2))'
	
	immath(imagename = [image[i], image[i+1]], mode='evalexpr',expr = expr, outfile = pre + "_spix_err_map_" +  str(int(freq[i])) + "_" + str(int(freq[i+1])) + "MHz_" + suf)
	
	exportfits(imagename = pre + "_spix_err_map_" +  str(int(freq[i])) + "_" + str(int(freq[i+1])) + "MHz_" + suf, fitsimage = pre + "_spix_err_map_" +  str(int(freq[i])) + "_" + str(int(freq[i+1])) + "MHz_" + suf + ".fits", overwrite=True)
	
	spix_err_im.append(pre + "_spix_err_map_" +  str(int(freq[i])) + "_" + str(int(freq[i+1])) + "MHz_" + suf)
	
#	shutil.rmtree(pre + "_spix_err_map_" +  str(int(freq[i])) + "_" + str(int(freq[i+1])) + "MHz_" + suf)

# Spectral Curvature Map
if len(freq) >= 3 :
	for i in range(len(spix_im)-1):
		expr = '-IM0 + IM1'
		immath(imagename = [spix_im[i], spix_im[i+1]], mode='evalexpr',expr = expr, outfile = pre + "_spix_curve_map_" + str(int(freq[i])) + "_" + str(int(freq[i+2])) + "MHz_" + suf)
		exportfits(imagename = pre + "_spix_curve_map_" + str(int(freq[i])) + "_" + str(int(freq[i+2])) + "MHz_" + suf, fitsimage = pre + "_spix_curve_map_" + str(int(freq[i])) + "_" + str(int(freq[i+2])) + "MHz_" + suf + ".fits", overwrite=True)
		shutil.rmtree(pre + "_spix_curve_map_" + str(int(freq[i])) + "_" + str(int(freq[i+2])) + "MHz_" + suf)
		
		expr = 'sqrt(pow(0.00009/IM0,2) + pow(-0.0000275/IM1,2) + pow(0.00000773/IM2,2))/0.46679'
		immath(imagename = [image[0], image[1], image[2]], mode='evalexpr',expr = expr, outfile = pre + "_spix_curve_err_map_" + str(int(freq[i])) + "_" + str(int(freq[i+2])) + "MHz_" + suf)
		exportfits(imagename = pre + "_spix_curve_err_map_" + str(int(freq[i])) + "_" + str(int(freq[i+2])) + "MHz_" + suf, fitsimage = pre + "_spix_curve_err_map_" + str(int(freq[i])) + "_" + str(int(freq[i+2])) + "MHz_" + suf + ".fits", overwrite=True)
		shutil.rmtree(pre + "_spix_curve_err_map_" + str(int(freq[i])) + "_" + str(int(freq[i+2])) + "MHz_" + suf)
		
	
# Remove extra directories
for i in image:
	if i == ref_image : 
		shutil.rmtree(ref_image)
	else :
		shutil.rmtree(i)
		shutil.rmtree(i.replace(".regrid", ""))

for i in spix_im:
	shutil.rmtree(i)
	
for i in spix_err_im:
	shutil.rmtree(i)

print(' ')
print('###########################   THE END   #########################')
print('Enjoy!, ... or Not???')
# End

