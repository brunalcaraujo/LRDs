import numpy as np
from astropy.io import fits
import warnings

def fnu_to_flambda(fnu, wave, wave_unit="um"):
    """
    Convert F_nu to F_lambda.

    Parameters
    ----------
    fnu : array_like
        Flux density per frequency [erg/s/cm^2/Hz]
    wave : array_like
        Wavelength
    wave_unit : {'um', 'A'}
        Unit of wavelength

    Returns
    -------
    flambda : ndarray
        Flux density per wavelength [erg/s/cm^2/Å]
    """

    c = 2.998e10  # cm/s

    wave = np.asarray(wave)

    if wave_unit == "um":
        lam_cm = wave * 1e-4
    elif wave_unit == "A":
        lam_cm = wave * 1e-8
    else:
        raise ValueError("wave_unit must be 'um' or 'A'")

    return fnu * c / lam_cm**2 / 1e8

def to_restframe(wave, flux, z, flux_type="flambda"):
    """
    Convert spectrum to rest frame.

    Parameters
    ----------
    wave : array_like
        Observed wavelength
    flux : array_like
        Observed flux density
    z : float
        Redshift
    flux_type : {'flambda', 'fnu'}

    Returns
    -------
    wave_rest, flux_rest
    """

    z1 = 1.0 + z
    wave_rest = wave / z1

    if flux_type == "flambda":
        flux_rest = flux * z1
    elif flux_type == "fnu":
        flux_rest = flux / z1
    else:
        raise ValueError("flux_type must be 'flambda' or 'fnu'")

    return wave_rest, flux_rest


def read_spectrum_fits(
    fits_path,
    wave_col="wave",
    flux_col="flux"
):
    """
    Read wavelength and flux from a FITS spectrum.

    Returns
    -------
    wave : ndarray
    flux : ndarray
    """

    with fits.open(fits_path) as hdul:
        data = hdul[1].data
        wave = np.asarray(data[wave_col])
        flux = np.asarray(data[flux_col])

    return wave, flux

def normalize_spectrum(
    wave,
    flux,
    window=(0.3546,0.3746), #ao redor de 3646
    statistic="median",
    min_points=4
):
    """
    Normalize a spectrum using a continuum window.

    Parameters
    ----------
    wave : array_like
        Wavelength array (rest-frame)
    flux : array_like
        Flux array
    window : tuple
        (lambda_min, lambda_max) in same units as wave
    statistic : {'median', 'mean'}
        Statistic used to estimate continuum level
    min_points : int
        Minimum number of points required in the window

    Returns
    -------
    flux_norm : ndarray
        Normalized flux
    norm_factor : float
        Normalization factor applied (original flux was divided by this)
    """

    wave = np.asarray(wave)
    flux = np.asarray(flux)

    mask = (wave >= window[0]) & (wave <= window[1])

    if mask.sum() < min_points:
        raise ValueError(
            f"Not enough points in normalization window {window}"
        )

    if statistic == "median":
        norm_factor = np.nanmedian(flux[mask])
    elif statistic == "mean":
        norm_factor = np.nanmean(flux[mask])
    else:
        raise ValueError("statistic must be 'median' or 'mean'")

    if not np.isfinite(norm_factor) or norm_factor <= 0:
        raise ValueError("Invalid normalization factor")

    flux_norm = flux / norm_factor

    return flux_norm, norm_factor


def load_spectrum(
    fits_path,
    z=None,
    input_flux_unit="uJy",
    wave_unit="um",
    restframe=True,
    normalize=False,
    norm_window=(0.3546,0.3746), #ao redor de 3646
    norm_statistic="median",
    output_flux_scale=None, 
):
    """
    Load a spectrum, convert units, optionally shift to rest frame
    and normalize.
    """

    wave, flux = read_spectrum_fits(fits_path)

    # ---- Flux unit conversion ----
    if input_flux_unit == "uJy":
        fnu = flux * 1e-29
        flambda = fnu_to_flambda(fnu, wave, wave_unit)
        flux = flambda
        flux_type = "flambda"
    else:
        raise NotImplementedError("Only uJy implemented for now")

    # ---- Rest-frame correction ----
    if restframe and z is not None:
        wave, flux = to_restframe(wave, flux, z, flux_type)

    # ---- Normalization ----
    norm_factor = None
    normalized = False
    norm_error = None

    if normalize:
        try:
            flux, norm_factor = normalize_spectrum(
                wave,
                flux,
                window=norm_window,
                statistic=norm_statistic,
            )
            normalized = True

        except ValueError as e:
            # Store error information, but DO NOT crash
            norm_error = str(e)
            normalized = False

    # ---- Optional scaling (for plotting convenience) ----
    if output_flux_scale is not None:
        flux = flux * output_flux_scale

    return {
        "wave": wave,
        "flux": flux,
        "z": z,
        "file": fits_path,
        "normalized": normalized,
        "norm_window": norm_window if normalized else None,
        "norm_factor": norm_factor,
        "norm_error": norm_error,
        "output_flux_scale": output_flux_scale,
    }