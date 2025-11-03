import numpy as np
from astropy.io import fits
import logging
import os

# Input FITS images
IMAGE_FILES = [
    "A85_reproj_A85_cutout_A85_400MHz_spix_100luvcut_rob-0.5-MFS-image_bm715631_pbcor.fits",
    "A85_reproj_A85_cutout_A85_700MHz_spix_100luvcut_rob-0.5-MFS-image_bm715631_pbcor.fits",
    "A85_reproj_A85_cutout_A85_1280MHz_spix_100luvcut_rob-0.5_nt2-MFS-image_bm715631_PBCOR.fits"
]

# Alignment file and manual correction
# Finetune Re-align image; [x, y] pixels, +x to shift right and +y to shift up
XY_CORR = [[1, -2], [1, -1]] 

#----------------------
def realign_images(images, xy_corr, output_prefix="realign"):
    """Shift images according to xy_corr and update FITS headers."""
    aligned = [images[0]]  # first image is reference, no shift

    for i in range(len(images) - 1):
        f = images[i + 1]
        dx, dy = map(int, np.round(xy_corr[i]))

        # Open image
        hdu = fits.open(f)[0]
        data = np.squeeze(hdu.data)  # remove singleton axes
        header = hdu.header

        # Determine dimensionality
        if data.ndim == 2:
            shifted = np.roll(np.roll(data, dy, axis=0), dx, axis=1)
        elif data.ndim == 3:
            shifted = np.roll(np.roll(data, dy, axis=1), dx, axis=2)
        elif data.ndim == 4:
            shifted = np.roll(np.roll(data, dy, axis=2), dx, axis=3)
        else:
            raise ValueError(f"Unsupported FITS data shape: {data.shape}")

        outname = f"{output_prefix}_{os.path.basename(f)}"
        fits.writeto(outname, shifted, header, overwrite=True)
        logging.info(f"Aligned: {outname}")
        aligned.append(outname)

    return aligned


# -------------------------------------------------------------------
# MAIN EXECUTION
# -------------------------------------------------------------------
if __name__ == "__main__":

    # Step 2: Realign images
    aligned_files = realign_images(IMAGE_FILES, XY_CORR)
    
    logging.info("Re-alignment completed successfully.")
    
