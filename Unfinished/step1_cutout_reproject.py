from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D
from reproject import reproject_interp
import astropy.units as u
import numpy as np
import os
from astropy.coordinates import Angle

def radec_to_degrees(ra_str, dec_str):
    """
    Convert RA/DEC from sexagesimal to decimal degrees.

    Parameters
    ----------
    ra_str : str
        RA in 'hh:mm:ss.s' format
    dec_str : str
        DEC in 'dd:mm:ss.s' format

    Returns
    -------
    tuple
        (RA_deg, DEC_deg) in decimal degrees
    """
    ra_deg = Angle(ra_str, unit='hourangle').deg
    dec_deg = Angle(dec_str, unit='deg').deg
    return ra_deg, dec_deg


def cutout_images_radec(files, output_prefix, ra, dec, size_arcmin):
    """
    Cutout multiple FITS images centered at (RA, Dec) with size in arcminutes.

    Parameters
    ----------
    files : list of str
        List of FITS files
    output_prefix : str
        Prefix for output files
    ra, dec : float
        RA, Dec in degrees
    size_arcmin : float or tuple
        Cutout size in arcminutes (float = square, tuple = (x, y))

    Returns
    -------
    cutout_files : list of str
        List of output cutout FITS files
    """
    cutout_files = []

    for f in files:
        # Open input FITS
        with fits.open(f) as hdul:
            data = hdul[0].data
            header = hdul[0].header

        # Collapse to 2D
        data = np.squeeze(data)
        if data.ndim > 2:
            data = data[0]

        # Extract WCS (celestial only)
        wcs = WCS(header).celestial

        # Define center position
        position = SkyCoord(ra, dec, unit="deg", frame="icrs")

        # Define cutout size
        if isinstance(size_arcmin, (list, tuple)):
            size_quantity = (size_arcmin[0] * u.arcmin, size_arcmin[1] * u.arcmin)
        else:
            size_quantity = (size_arcmin * u.arcmin, size_arcmin * u.arcmin)

        # Make cutout (allow partial overlap, fill with NaN if outside)
        cutout = Cutout2D(
            data=data,
            position=position,
            size=size_quantity,
            wcs=wcs,
            mode="partial",
            fill_value=np.nan
        )

        # New header from cutout WCS
        header_cut = cutout.wcs.to_header()
        header_cut["NAXIS"]  = 2
        header_cut["NAXIS1"] = cutout.data.shape[1]
        header_cut["NAXIS2"] = cutout.data.shape[0]

        # Preserve some metadata
        for key in ["BUNIT", "BMAJ", "BMIN", "BPA"]:
            if key in header:
                header_cut[key] = header[key]

        # Output file name
        out_file = f"{output_prefix}_cutout_{os.path.basename(f)}"

        # Write cutout
        fits.writeto(out_file, cutout.data, header_cut, overwrite=True)

        print(f"Saved cutout: {out_file}")
        cutout_files.append(out_file)

    return cutout_files


def reproject_images(files, reference_index=0, output_prefix="reproj"):
    """
    Reproject multiple FITS images onto the grid of a reference image.

    Parameters
    ----------
    files : list of str
        List of FITS files (already cutouts)
    reference_index : int
        Index of file to use as reference grid
    output_prefix : str
        Prefix for output files

    Returns
    -------
    reprojected_files : list of str
        List of output reprojected FITS files
    """
    reprojected_files = []

    # Reference image
    ref_file = files[reference_index]
    with fits.open(ref_file) as hdul:
        ref_data = np.squeeze(hdul[0].data)
        ref_header = hdul[0].header
        ref_wcs = WCS(ref_header).celestial

    # Clean reference header
    header_out = ref_wcs.to_header()
    header_out["NAXIS"]  = 2
    header_out["NAXIS1"] = ref_data.shape[-1]
    header_out["NAXIS2"] = ref_data.shape[-2]

    for i, f in enumerate(files):
        with fits.open(f) as hdul:
            data = np.squeeze(hdul[0].data)
            wcs = WCS(hdul[0].header).celestial
            hdu_clean = fits.PrimaryHDU(data=data, header=wcs.to_header())

        if i == reference_index:
            # Just copy reference
            out_file = f"{output_prefix}_{os.path.basename(f)}"
            fits.writeto(out_file, ref_data, ref_header, overwrite=True)
            print(f"Reference image saved as: {out_file}")
            reprojected_files.append(out_file)
            continue

        # Reproject
        array, footprint = reproject_interp(
            hdu_clean, header_out,
            shape_out=(ref_data.shape[-2], ref_data.shape[-1])
        )

        # Save with reference header
        out_file = f"{output_prefix}_{os.path.basename(f)}"
        fits.writeto(out_file, array, ref_header, overwrite=True)

        print(f"Reprojected image saved: {out_file}")
        reprojected_files.append(out_file)

    return reprojected_files


# ========================
# Example usage
# ========================
if __name__ == "__main__":
    IMAGE_FILES = [
        "A85_400MHz_spix_100luvcut_rob-0.5-MFS-image_bm715631_pbcor.fits",
        "A85_700MHz_spix_100luvcut_rob-0.5-MFS-image_bm715631_pbcor.fits",
        "A85_1280MHz_spix_100luvcut_rob-0.5_nt2-MFS-image_bm715631_PBCOR.fits"
    ]

#    # Input position (degrees)
#    RA_DEG = 10.445   # example
#    DEC_DEG = -9.183  # example
#    SIZE_ARCMIN = 10  # cutout size
    
    # Cropping options
    # RA/DEC in sexagesimal
    ra_hms = "00:41:51.3"
    dec_dms = "-09:19:09"
    SIZE_ARCMIN = (10.0, 10.0)         # 5 arcmin x 5 arcmin
    
    # Convert to degrees
    center_deg = radec_to_degrees(ra_hms, dec_dms)
    print("Center (deg):", center_deg)
    RA_DEG, DEC_DEG = center_deg

    # Step 1: Cutouts
    cutout_files = cutout_images_radec(
        IMAGE_FILES,
        output_prefix="A85",
        ra=RA_DEG,
        dec=DEC_DEG,
        size_arcmin=SIZE_ARCMIN
    )

    # Step 2: Reproject cutouts
    reprojected_files = reproject_images(
        cutout_files,
        reference_index=0,   # first image is reference
        output_prefix="A85_reproj"
    )

