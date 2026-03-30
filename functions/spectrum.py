import numpy as np
from astropy.io import fits
from scipy.signal import savgol_filter
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
    flux_col="flux",
    err_col="err"
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
        err = np.asarray(data[err_col])

    return wave, flux, err

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


def smooth_spectrum(
    flux,
    method="savgol",
    window=11,
    polyorder=2,
):
    """
    Smooth a spectrum flux array.

    Parameters
    ----------
    flux : array_like
        Flux array
    method : {'savgol', 'boxcar'}
        Smoothing method
    window : int
        Window size (number of points)
    polyorder : int
        Polynomial order for Savitzky-Golay filter

    Returns
    -------
    flux_smooth : ndarray
        Smoothed flux
    """

    flux = np.asarray(flux)

    # Remove NaNs via interpolation
    if np.any(~np.isfinite(flux)):
        x = np.arange(len(flux))
        mask = np.isfinite(flux)
        if mask.sum() < 5:
            raise ValueError("Too many NaNs for smoothing")
        flux = np.interp(x, x[mask], flux[mask])

    n = len(flux)

    if window >= n:
        raise ValueError("Smoothing window larger than spectrum")

    if method == "savgol":
        if window % 2 == 0:
            window += 1
        return savgol_filter(flux, window_length=window, polyorder=polyorder)

    elif method == "boxcar":
        kernel = np.ones(window) / window
        return np.convolve(flux, kernel, mode="same")

    else:
        raise ValueError("method must be 'savgol' or 'boxcar'")

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
    smooth=False,
    smooth_method="savgol",
    smooth_window=11,
    smooth_polyorder=2,
 
):
    """
    Load a spectrum, convert units, optionally shift to rest frame
    and normalize.
    """

    wave, flux, err = read_spectrum_fits(fits_path)

    # ---- Flux unit conversion ----
    if input_flux_unit == "uJy":
        fnu = flux * 1e-29
        err_fnu = err * 1e-29

        flambda = fnu_to_flambda(fnu, wave, wave_unit)
        err_flambda = fnu_to_flambda(err_fnu, wave, wave_unit)

        flux = flambda
        err = err_flambda
        flux_type = "flambda"
    else:
        raise NotImplementedError("Only uJy implemented for now")

    # ---- Rest-frame correction ----
    if restframe and z is not None:
        wave, flux = to_restframe(wave, flux, z, flux_type)

        # aplicar MESMA transformação no erro
        z1 = 1.0 + z
        if flux_type == "flambda":
            err = err * z1
        elif flux_type == "fnu":
            err = err / z1

    # ---- Smoothing ----
    smoothed = False
    smooth_error = None

    if smooth:
        try:
            flux = smooth_spectrum(
                flux,
                method=smooth_method,
                window=smooth_window,
                polyorder=smooth_polyorder,
            )
            smoothed = True
        except ValueError as e:
            smooth_error = str(e)
            smoothed = False


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
            err = err / norm_factor #normalizando o erro também
            normalized = True

        except ValueError as e:
            # Store error information, but DO NOT crash
            norm_error = str(e)
            normalized = False

    # ---- Optional scaling (for plotting convenience) ----
    if output_flux_scale is not None:
        flux = flux * output_flux_scale

    # ---- Noise calculation ----
    noise_metrics = compute_noise_metrics(wave, flux, err)

    return {
        "wave": wave,
        "flux": flux,
        "z": z,
        "file": fits_path,
        "normalized": normalized,
        "norm_window": norm_window if normalized else None,
        "norm_factor": norm_factor,
        "norm_error": norm_error,
        "smoothed": smoothed,
        "smooth_method": smooth_method if smoothed else None,
        "smooth_window": smooth_window if smoothed else None,
        "smooth_error": smooth_error,
        "output_flux_scale": output_flux_scale,
        "noise": noise_metrics,
    }


def compute_noise_metrics(
    wave,
    flux,
    err=None,
    smooth=True,
    smooth_window=21,
    smooth_polyorder=2,
):
    """
    Compute noise metrics for a spectrum.

    Parameters
    ----------
    wave : array
    flux : array
    err : array or None
        Per-pixel uncertainty
    smooth : bool
        Whether to estimate empirical noise
    """

    wave = np.asarray(wave)
    flux = np.asarray(flux)

    # Remove invalid values
    mask = np.isfinite(flux)
    if err is not None:
        err = np.asarray(err)
        mask &= np.isfinite(err) & (err > 0)

    flux = flux[mask]
    if err is not None:
        err = err[mask]

    results = {}

    # -----------------------------
    # (A) Noise from error array
    # -----------------------------
    if err is not None:
        noise_median = np.median(err)
        snr = flux / err
        snr_median = np.median(snr)
    else:
        noise_median = np.nan
        snr_median = np.nan

    results["noise_from_err"] = noise_median
    results["snr_median"] = snr_median

    # -----------------------------
    # (B) Empirical noise
    # -----------------------------
    if smooth:
        try:
            flux_smooth = smooth_spectrum(
                flux,
                method="savgol",
                window=smooth_window,
                polyorder=smooth_polyorder,
            )

            residuals = flux - flux_smooth
            noise_empirical = np.std(residuals)

        except Exception:
            noise_empirical = np.nan
    else:
        noise_empirical = np.nan

    results["noise_empirical"] = noise_empirical

    # -----------------------------
    # (C) Robust SNR (median-based)
    # -----------------------------
    if err is not None:
        snr_robust = np.nanmedian(flux) / np.nanmedian(err)
    else:
        snr_robust = np.nan

    results["snr_robust"] = snr_robust


    # -----------------------------
    # (D) Empirical SNR
    # -----------------------------
    if np.isfinite(noise_empirical) and noise_empirical > 0:
        snr_empirical = np.nanmedian(flux) / noise_empirical
    else:
        snr_empirical = np.nan

    results["snr_empirical"] = snr_empirical

    return results