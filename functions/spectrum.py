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
    err : ndarray
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
    err=None,
    window=(0.3446,0.3646), #antes de 3646
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
    err : array_like
        Error array
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

    if err is not None:
        err = np.asarray(err)

    mask = (wave >= window[0]) & (wave <= window[1])

    if mask.sum() < min_points:
        raise ValueError(
            f"Not enough points in normalization window {window}"
        )

    # ---- continuum level -----
    if statistic == "median":
        norm_factor = np.nanmedian(flux[mask])
    elif statistic == "mean":
        norm_factor = np.nanmean(flux[mask])
    else:
        raise ValueError("statistic must be 'median' or 'mean'")

    if not np.isfinite(norm_factor) or norm_factor <= 0:
        raise ValueError("Invalid normalization factor")

    flux_norm = flux / norm_factor

    # ---- normalize error ----
    if err is not None:
        err_norm = err / norm_factor

        # estatística do erro na janela
        err_window_mean = np.nanmean(err_norm[mask])
        err_window_median = np.nanmedian(err_norm[mask])
    else:
        err_norm = None
        err_window_mean = None
        err_window_median = None

    return flux_norm, err_norm, norm_factor, err_window_mean, err_window_median



def load_spectrum(
    fits_path,
    z=None,
    input_flux_unit="uJy",
    wave_unit="um",
    restframe=True,
    normalize=False,
    norm_window=(0.3446,0.3646), #antes de 3646
    norm_statistic="median",
    output_flux_scale=None,
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


    # ---- Normalization ----
    norm_factor = None
    normalized = False
    norm_error = None

    if normalize:
        try:
            flux, err, norm_factor, err_mean, err_median = normalize_spectrum(
                wave,
                flux,
                err=err,
                window=norm_window,
                statistic=norm_statistic,
            )

            normalized = True

        except ValueError as e:
            norm_error = str(e)
            normalized = False
            err_mean = None
            err_median = None

    # ---- Optional scaling (for plotting convenience) ----
    if output_flux_scale is not None:
        flux = flux * output_flux_scale

    return {
        "wave": wave,
        "flux": flux,
        "err": err,
        "z": z,
        "file": fits_path,
        "normalized": normalized,
        "norm_window": norm_window if normalized else None,
        "norm_factor": norm_factor,
        "norm_error": norm_error,
        "norm_err_mean": err_mean if normalized else None,
        "norm_err_median": err_median if normalized else None,
        "output_flux_scale": output_flux_scale,
    }

def compute_error_stats(
    wave,
    flux,
    err,
    window,
    normalized=True,   # 🔥 controla o comportamento
    statistic="both",
    min_points=3
):
    """
    Compute flux and error statistics in a given wavelength window.

    Parameters
    ----------
    wave : array_like
        Wavelength array
    flux : array_like
        Flux array (preferably normalized)
    err : array_like
        Error array (same units as flux)
    window : tuple
        (lambda_min, lambda_max)
    statistic : {'mean', 'median', 'both'}
        Which statistics to return
    min_points : int
        Minimum number of points required

    Returns
    -------
    dict with statistics
    """

    wave = np.asarray(wave)
    flux = np.asarray(flux)
    err = np.asarray(err)

    mask = (wave >= window[0]) & (wave <= window[1])
    n = mask.sum()

    # Caso NÃO esteja normalizado → tudo vira NaN
    if not normalized:
        return {
            "n_points": n,
            "flux_mean": np.nan,
            "flux_median": np.nan,
            "err_mean": np.nan,
            "err_median": np.nan,
            "err_rms": np.nan,
            "snr": np.nan
        }

    # poucos pontos
    if n < min_points:
        return {
            "n_points": n,
            "flux_mean": np.nan,
            "flux_median": np.nan,
            "err_mean": np.nan,
            "err_median": np.nan,
            "err_rms": np.nan,
            "snr": np.nan
        }

    flux_window = flux[mask]
    err_window = err[mask]

    results = {
        "n_points": n,
    }

    # ---- Flux stats ----
    if statistic in ["mean", "both"]:
        results["flux_mean"] = np.nanmean(flux_window)

    if statistic in ["median", "both"]:
        results["flux_median"] = np.nanmedian(flux_window)

    # ---- Error stats ----
    if statistic in ["mean", "both"]:
        results["err_mean"] = np.nanmean(err_window)

    if statistic in ["median", "both"]:
        results["err_median"] = np.nanmedian(err_window)

    # ---- RMS error (mais físico) ----
    err_rms = np.sqrt(np.nanmean(err_window**2))
    results["err_rms"] = err_rms

    # ---- SNR (definição simples e útil) ----
    # usa fluxo médio / erro rms
    if "flux_mean" in results and err_rms > 0:
        results["snr"] = results["flux_mean"] / err_rms
    else:
        results["snr"] = np.nan

    return results

def compute_mean_spectrum(
    spectra_list,
    wave_grid=None,
    n_clip_end=0,
    interp_kind="linear",
    return_error=True,
    flux_min=None
):
    """
    Compute mean spectrum from a list of spectra.

    Parameters
    ----------
    spectra_list : list of dict
        Output of load_spectrum()
    wave_grid : array_like or None
        Common wavelength grid. If None, will be auto-generated.
    n_clip_end : int
        Number of points to remove from the END of each spectrum
    interp_kind : str
        Interpolation type (currently only linear implemented)
        ---- Later I can add other type if necessary 
    return_error : bool
        If True, compute error on the mean

    Returns
    -------
    dict with:
        wave
        flux_mean
        flux_std
        err_mean (optional)
        n_contrib (number of spectra contributing per pixel)
        n_objects (number of spectra used)
    """

    # -------------------------
    # 1. construir grade comum
    # -------------------------
    if wave_grid is None:
        # máximos e mínimos de cada espectro
        wmins = [np.nanmin(s["wave"]) for s in spectra_list]
        wmaxs = [np.nanmax(s["wave"]) for s in spectra_list]

        # cobertura total
        wmin = np.nanmin(wmins)
        wmax = np.nanmax(wmaxs)

        # resolução baseada no primeiro espectro
        ref_wave = spectra_list[0]["wave"]
        dw = np.nanmedian(np.diff(ref_wave))

        wave_grid = np.arange(wmin, wmax, dw)

    wave_grid = np.asarray(wave_grid)

    # -------------------------
    # 2. interpolar espectros
    # -------------------------
    flux_stack = []
    err_stack = []

    for spec in spectra_list:
        wave = spec["wave"]
        flux = spec["flux"]
        err = spec["err"]

        # cortar final ruidoso
        if n_clip_end > 0:
            wave = wave[:-n_clip_end]
            flux = flux[:-n_clip_end]
            if err is not None:
                err = err[:-n_clip_end]

        # interpolação
        flux_interp = np.interp(
            wave_grid,
            wave,
            flux,
            left=np.nan,
            right=np.nan
        )

        flux_stack.append(flux_interp)

        if return_error and err is not None:
            err_interp = np.interp(
                wave_grid,
                wave,
                err,
                left=np.nan,
                right=np.nan
            )
            err_stack.append(err_interp)

    flux_stack = np.array(flux_stack)

    # -------------------------
    # limpeza de valores ruins
    # -------------------------
    if flux_min is not None:
        flux_stack[flux_stack < flux_min] = np.nan

    n_objects = flux_stack.shape[0]

    # -------------------------
    # 3. média ignorando NaN
    # -------------------------
    flux_mean = np.nanmean(flux_stack, axis=0)
    flux_std = np.nanstd(flux_stack, axis=0)

    # quantos espectros contribuíram em cada pixel
    n_contrib = np.sum(~np.isnan(flux_stack), axis=0)

    results = {
        "wave": wave_grid,
        "flux_mean": flux_mean,
        "flux_std": flux_std,
        "n_contrib": n_contrib,
        "n_objects": n_objects 
    }

    # -------------------------
    # 4. erro do espectro médio
    # -------------------------
    if return_error and len(err_stack) > 0:
        err_stack = np.array(err_stack)

        if flux_min is not None:
            err_stack[flux_stack < flux_min] = np.nan

        # erro médio propagado (aproximação)
        err_mean = np.sqrt(np.nanmean(err_stack**2, axis=0)) / np.sqrt(n_contrib)

        results["err_mean"] = err_mean

    return results
