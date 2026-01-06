import numpy as np
from astropy.io import fits

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

def load_spectrum(
    fits_path,
    z=None,
    input_flux_unit="uJy",
    wave_unit="um",
    restframe=True,
    output_flux_scale=None,
):
    """
    Load a spectrum, convert units, and optionally shift to rest frame.

    Returns
    -------
    spectrum : dict
        Keys: wave, flux, z
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

    # ---- Optional scaling (for plotting convenience) ----
    if output_flux_scale is not None:
        flux = flux * output_flux_scale

    return {
        "wave": wave,
        "flux": flux,
        "z": z,
        "file": fits_path
    }