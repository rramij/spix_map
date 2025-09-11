#!/usr/bin/env python3
"""
Spectral Index and Curvature Map Generator: Upto 3 frequencies
Author: Ramij Raja
Date: 11 Sep 2025

Dependencies:
    numpy >= 1.12.0
    astropy
"""

import os
import shutil
import logging
import numpy as np
from astropy.io import fits

# -------------------------------------------------------------------
# USER INPUTS (edit these as needed)
# -------------------------------------------------------------------

# Input FITS images
IMAGE_FILES = [
    "A85_400MHz_spix_100luvcut_rob-0.5-MFS-image_bm715631_pbcor.fits",
    "A85_700MHz_spix_100luvcut_rob-0.5-MFS-image_bm715631_pbcor.fits",
    "A85_1280MHz_spix_100luvcut_rob-0.5_nt2-MFS-image_bm715631_PBCOR.fits"
]

# Alignment file and manual correction
ALIGN_FILE = "align.dat"
EXTRA_CORR = [[0, 0], [1, 0]]   # [dx, dy] in pixels, applied after align.dat

# Observing frequencies (MHz)
FREQS = [400., 700., 1280.]

# Image RMS noise levels (Jy)
RMS = [5.0e-5, 2e-5, 7e-6]

# Fractional flux density errors
F_ERR = [0.1, 0.07, 0.05]

# Masking threshold: N × rms
SIGMA_LEVEL = 3

# Reference frequency for headers (keep same as first freq)
REF_FREQ = FREQS[0]

# Output filenames
OUTPUT_PREFIX = "A85"   # e.g., Cluster name
OUTPUT_SUFFIX = "B7x5"  # e.g., beam resolution

# Central cropped image size (pixels)
IMSIZE = 2000

# Optional spectral index masking range
MASK_RANGE = False       # True/False
MASK_PIX = [-5.0, 0.0]   # [min, max]

# Keep intermediate files?
KEEP_INTERMEDIATE = False   # if True, move to "intermediate/" folder

# -------------------------------------------------------------------
# Logging setup
# -------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s"
)

# -------------------------------------------------------------------
# Utility functions
# -------------------------------------------------------------------
def read_alignment(file, extra_corr=None):
    """Read alignment offsets from file and apply extra correction."""
    data = np.loadtxt(file, skiprows=1)
    ref = data[0]
    offsets = ref[1:] - data[1:, 1:]
    if extra_corr is not None:
        offsets += np.array(extra_corr)
    return offsets


def realign_images(images, xy_corr):
    """Shift images according to xy_corr."""
    aligned = [images[0]]  # first one is reference
    for i in range(len(images) - 1):
        data, hdr = fits.getdata(images[i + 1], header=True)
        dx, dy = map(int, np.round(xy_corr[i]))
        shifted = np.roll(data, (dx, dy), axis=(3, 2))
        outname = images[i + 1].replace(".fits", "_realign.fits")
        fits.writeto(outname, shifted, hdr, overwrite=True)
        logging.info(f"Aligned: {outname}")
        aligned.append(outname)
    return aligned


def mask_and_crop(image_file, rms, sig, imsize):
    """Mask pixels below sigma*rms and crop to central region."""
    data, hdr = fits.getdata(image_file, header=True)
    mask_level = sig * rms
    data[data < mask_level] = np.nan

    x_center, y_center = int(hdr["CRPIX1"]), int(hdr["CRPIX2"])
    x1, x2 = x_center - imsize // 2, x_center + imsize // 2
    y1, y2 = y_center - imsize // 2, y_center + imsize // 2

    cropped = np.full_like(data, np.nan)
    cropped[..., y1:y2, x1:x2] = data[..., y1:y2, x1:x2]

    outname = image_file.replace(".fits", "_masked.fits")
    fits.writeto(outname, cropped, hdr, overwrite=True)
    return outname


def spectral_index(data0, data1, f0, f1, rms0, rms1, ferr0, ferr1):
    """Compute spectral index and its error."""
    alpha = np.log(data0 / data1) / np.log(f0 / f1)

    # Optional masking range
    if MASK_RANGE:
        alpha[(alpha <= min(MASK_PIX)) | (alpha >= max(MASK_PIX))] = np.nan

    delta_S0 = np.sqrt(rms0**2 + (ferr0 * data0) ** 2)
    delta_S1 = np.sqrt(rms1**2 + (ferr1 * data1) ** 2)
    delta_alpha = (1.0 / abs(np.log(f0 / f1))) * np.sqrt(
        (delta_S0 / data0) ** 2 + (delta_S1 / data1) ** 2
    )

    mask = np.isnan(alpha)
    delta_alpha[mask] = np.nan
    return alpha, delta_alpha


def curvature(spix0, spix1, freqs, rms, images):
    """Compute curvature and error map."""
    curv = spix1 - spix0
    f0, f1, f2 = freqs
    im0, im1, im2 = images
    a = -np.log(f1 / f2) * rms[0]
    b = np.log(f0 / f2) * rms[1]
    c = -np.log(f0 / f1) * rms[2]
    d = np.log(f0 / f1) * np.log(f1 / f2)
    curv_e = np.sqrt((a / im0) ** 2 + (b / im1) ** 2 + (c / im2) ** 2) / d
    curv_e[np.isnan(curv)] = np.nan
    return curv, curv_e


def cleanup(files, extensions=("_realign.fits", "_masked.fits", "_realign_masked.fits")):
    """Remove or archive intermediate files."""
    if KEEP_INTERMEDIATE:
        folder = "intermediate"
        os.makedirs(folder, exist_ok=True)
        for f in files:
            for ext in extensions:
                candidate = f.replace(".fits", ext)
                if os.path.exists(candidate):
                    shutil.move(candidate, os.path.join(folder, os.path.basename(candidate)))
                    logging.info(f"Moved intermediate to {folder}/: {candidate}")
    else:
        for f in files:
            for ext in extensions:
                candidate = f.replace(".fits", ext)
                if os.path.exists(candidate):
                    os.remove(candidate)
                    logging.info(f"Removed intermediate: {candidate}")


# -------------------------------------------------------------------
# Main pipeline
# -------------------------------------------------------------------
def main():
    logging.info("Starting pipeline...")

    # Step 1: Alignment
    xy_corr = read_alignment(ALIGN_FILE, extra_corr=EXTRA_CORR)
    aligned = realign_images(IMAGE_FILES, xy_corr)

    # Step 2: Mask and crop
    masked = [
        mask_and_crop(aligned[i], RMS[i], SIGMA_LEVEL, IMSIZE)
        for i in range(len(aligned))
    ]

    # Step 3: Spectral Index Maps
    ref_hdr = fits.getheader(masked[0])
    spix_files, spix_err_files = [], []

    for i in range(len(masked) - 1):
        d0 = fits.getdata(masked[i])
        d1 = fits.getdata(masked[i + 1])
        alpha, alpha_err = spectral_index(
            d0, d1, FREQS[i], FREQS[i + 1], RMS[i], RMS[i + 1], F_ERR[i], F_ERR[i + 1]
        )

        out_spix = f"{OUTPUT_PREFIX}_spix_map_{int(FREQS[i])}_{int(FREQS[i+1])}MHz_{OUTPUT_SUFFIX}.fits"
        out_err = f"{OUTPUT_PREFIX}_spix_err_map_{int(FREQS[i])}_{int(FREQS[i+1])}MHz_{OUTPUT_SUFFIX}.fits"
        fits.writeto(out_spix, alpha, ref_hdr, overwrite=True)
        fits.writeto(out_err, alpha_err, ref_hdr, overwrite=True)

        logging.info(f"Saved: {out_spix}, {out_err}")
        spix_files.append(out_spix)
        spix_err_files.append(out_err)

    # Step 4: Curvature Maps
    for i in range(len(spix_files) - 1):
        s0 = fits.getdata(spix_files[i])
        s1 = fits.getdata(spix_files[i + 1])
        im0 = fits.getdata(masked[i])
        im1 = fits.getdata(masked[i + 1])
        im2 = fits.getdata(masked[i + 2])
        curv, curv_e = curvature(
            s0, s1, FREQS[i : i + 3], RMS[i : i + 3], [im0, im1, im2]
        )

        out_curv = f"{OUTPUT_PREFIX}_curve_map_{int(FREQS[i])}_{int(FREQS[i+2])}MHz_{OUTPUT_SUFFIX}.fits"
        out_curv_err = f"{OUTPUT_PREFIX}_curve_err_map_{int(FREQS[i])}_{int(FREQS[i+2])}MHz_{OUTPUT_SUFFIX}.fits"
        fits.writeto(out_curv, curv, ref_hdr, overwrite=True)
        fits.writeto(out_curv_err, curv_e, ref_hdr, overwrite=True)
        logging.info(f"Saved: {out_curv}, {out_curv_err}")

    # Step 5: Cleanup / Archive
    cleanup(IMAGE_FILES)
    logging.info("Pipeline completed successfully.")


if __name__ == "__main__":
    main()

