import numpy as np
from astropy.io import fits
import logging
import os
import shutil
from scipy.optimize import curve_fit
import multiprocessing as mp

# =====================================================
# ------------ USER INPUT PARAMETERS ------------------
# =====================================================
# Input FITS images
IMAGE_FILES = [
    "A85_reproj_A85_cutout_A85_400MHz_spix_100luvcut_rob-0.5-MFS-image_bm715631_pbcor.fits",
    "realign_A85_reproj_A85_cutout_A85_700MHz_spix_100luvcut_rob-0.5-MFS-image_bm715631_pbcor.fits",
    "realign_A85_reproj_A85_cutout_A85_1280MHz_spix_100luvcut_rob-0.5_nt2-MFS-image_bm715631_PBCOR.fits"
]


MASK_RANGE   = False                  # apply Spix masking?
MASK_PIX     = [-5.0, 0.0]            # valid spectral index range
KEEP_INTERMEDIATE = False             # keep intermediate files?
INTERMEDIATE_DIR  = "intermediate"    # folder for intermediates if kept
OUTPUT_PREFIX = "A85"              # prefix for output maps
OUTPUT_SUFFIX = "B7x5"                # suffix for output maps
NPROC         = 10                  # number of parallel processes (None=all cores)

# Frequencies in MHz (must match file order in align.dat)
FREQS    = np.array([400, 700, 1280])
RMS_LIST = [5.0e-5, 2e-5, 7e-6]  # rms per frequency (Jy/beam)
# Masking threshold: N × rms
SIGMA_LEVEL = 3
spix_freq = 500 # Calculate spectral index at frequency (MHz)

# =====================================================
# ------------ FUNCTIONS ------------------------------
# =====================================================

def spec_model(log_nu, A, B, C):
    """log S = A + B log ν + C (log ν)^2"""
    return A + B * log_nu + C * log_nu**2

def fit_pixel_weighted(args):
    """Fit one pixel with curve_fit using weights (1/rms), return alpha, beta, alpha_err, beta_err"""
    log_nu, fluxes, rms = args
    if np.any(~np.isfinite(fluxes)) or np.any(fluxes <= 0):
        return np.nan, np.nan, np.nan, np.nan

    log_S = np.log(fluxes)
    sigma = np.array(rms) / np.array(fluxes)  # σ_log = RMS / flux

    try:
        popt, pcov = curve_fit(spec_model, log_nu, log_S, sigma=sigma, absolute_sigma=True)
        _, B, C = popt
        # calculate spix
        spix = B + (2*C*np.log(spix_freq))
        perr = np.sqrt(np.diag(pcov))
        # calculate spix err
        spix_err = np.sqrt(perr[1]**2 + (2*np.log(spix_freq)*perr[2])**2)
        
        return spix, C, spix_err, perr[2]
    except Exception:
        return np.nan, np.nan, np.nan, np.nan

def fit_spectral_maps_curvefit_weighted(masked_files, freqs, rms_list,
                                        output_prefix, output_suffix,
                                        MASK_RANGE=False, MASK_PIX=[-5.0, 0.0],
                                        nproc=None,
                                        keep_intermediate=False,
                                        intermediate_dir="intermediate"):
    """Fit quadratic log-log spectrum per pixel with curve_fit + RMS weights + multiprocessing."""

    logging.info("Starting weighted spectral fitting with curve_fit + multiprocessing...")

    # --- Load data cube ---
    cube = []
    for f in masked_files:
        data = fits.getdata(f)
        data = np.squeeze(data)  # remove singleton axes
        if data.ndim != 2:
            raise ValueError(f"Image {f} is not 2D after squeezing (got shape {data.shape})")
        cube.append(data)

    cube = np.array(cube)  # shape: (n_freq, ny, nx)
    n_freq, ny, nx = cube.shape
    ref_hdr = fits.getheader(masked_files[0])

    # Flatten spatial grid
    cube = cube.reshape(n_freq, -1)  # (n_freq, n_pix)
    n_pix = cube.shape[1]

    log_nu = np.log(freqs)

    # --- Prepare args for parallel fitting ---
    args_list = [(log_nu, cube[:, i], rms_list) for i in range(n_pix)]

    # --- Parallel curve fitting per pixel ---
    with mp.Pool(processes=nproc) as pool:
        results = pool.map(fit_pixel_weighted, args_list)

    results = np.array(results)  # shape (n_pix, 4)

    # Unpack
    alpha     = results[:, 0].reshape(ny, nx)
    beta      = results[:, 1].reshape(ny, nx)
    alpha_err = results[:, 2].reshape(ny, nx)
    beta_err  = results[:, 3].reshape(ny, nx)

    # --- Apply spectral index masking ---
    if MASK_RANGE:
        logging.info(f"Applying spectral index mask: {MASK_PIX}")
        valid = (alpha >= MASK_PIX[0]) & (alpha <= MASK_PIX[1])
        alpha[~valid] = np.nan
        beta[~valid] = np.nan
        alpha_err[~valid] = np.nan
        beta_err[~valid] = np.nan

    # --- Save results ---
    out_spix     = f"{output_prefix}_spix_map_{int(freqs[0])}_{int(freqs[-1])}MHz_{output_suffix}.fits"
    out_spix_err = f"{output_prefix}_spix_err_map_{int(freqs[0])}_{int(freqs[-1])}MHz_{output_suffix}.fits"
    out_curv     = f"{output_prefix}_curve_map_{int(freqs[0])}_{int(freqs[-1])}MHz_{output_suffix}.fits"
    out_curv_err = f"{output_prefix}_curve_err_map_{int(freqs[0])}_{int(freqs[-1])}MHz_{output_suffix}.fits"

    fits.writeto(out_spix, alpha, ref_hdr, overwrite=True)
    fits.writeto(out_spix_err, alpha_err, ref_hdr, overwrite=True)
    fits.writeto(out_curv, beta, ref_hdr, overwrite=True)
    fits.writeto(out_curv_err, beta_err, ref_hdr, overwrite=True)

    logging.info(f"Saved spectral index maps: {out_spix}, {out_spix_err}")
    logging.info(f"Saved curvature maps: {out_curv}, {out_curv_err}")

    return out_spix, out_spix_err, out_curv, out_curv_err

def mask_sig(image_file, rms, sig):
    """Mask pixels below sigma*rms and crop to central region."""
    data, hdr = fits.getdata(image_file, header=True)
    mask_level = sig * rms
    data[data < mask_level] = np.nan

    outname = image_file.replace(".fits", "_masked.fits")
    fits.writeto(outname, data, hdr, overwrite=True)
    return outname
    
# =====================================================
# ------------------- MAIN PIPELINE -------------------
# =====================================================
if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    # Step 3: Use input FITS files
    raw_files = IMAGE_FILES
    
    # Step 2: Mask and crop
    masked_files = [
        mask_sig(raw_files[i], RMS_LIST[i], SIGMA_LEVEL)
        for i in range(len(raw_files))
    ]

    # Step 6: Spectral fitting
    fit_spectral_maps_curvefit_weighted(
        masked_files, FREQS, RMS_LIST,
        OUTPUT_PREFIX, OUTPUT_SUFFIX,
        MASK_RANGE=MASK_RANGE,
        MASK_PIX=MASK_PIX,
        nproc=NPROC,
        keep_intermediate=KEEP_INTERMEDIATE,
        intermediate_dir=INTERMEDIATE_DIR
    )

