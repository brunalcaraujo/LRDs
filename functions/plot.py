import matplotlib.pyplot as plt
import numpy as np
from .spectrum import load_spectrum


def plot_spectrum_ax(
    ax,
    spectrum,
    xlim=None,
    ylim=None,
    title=None,
    z=None,
    **plot_kwargs
):
    """
    Plot a spectrum on a given axis.

    Parameters
    ----------
    ax : matplotlib axis
    spectrum : dict
        Output of load_spectrum()
    xlim, ylim : tuple or None
        Axis limits, e.g. (xmin, xmax)
    title : str or None
    z : float or None
    plot_kwargs :
        Passed directly to ax.plot()
    """

    ax.plot(
        spectrum["wave"],
        spectrum["flux"],
        **plot_kwargs
    )

    if xlim is not None:
        ax.set_xlim(xlim)

    if ylim is not None:
        ax.set_ylim(ylim)

    ax.grid(True)

    if title is not None:
        ax.set_title(title, fontsize=10)

    if z is not None:
        ax.text(
            0.97, 0.95,
            f"z = {z:.3f}",
            transform=ax.transAxes,
            ha='right', va='top',
            fontsize=9,
            bbox=dict(facecolor='white', alpha=0.7, edgecolor='none')
        )

def make_spectrum_panel(
    spec_info,
    start=0,
    nrows=4,
    ncols=2,
    base_path="DeGraaff_espectros",
    xlim=None,
    ylim=None,
    plot_kwargs=None,
    loader_kwargs=None,
):
    """
    Parameters
    ----------
    spec_info : list of (filename, z)
    plot_kwargs : dict
        Passed to ax.plot()
    loader_kwargs : dict
        Passed to load_spectrum()
    """

    if plot_kwargs is None:
        plot_kwargs = {}

    if loader_kwargs is None:
        loader_kwargs = {}

    output_flux_scale = None

    n_panel = nrows * ncols
    subset = spec_info[start:start + n_panel]

    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(14, 4 * nrows),
        sharex=True, sharey=True
    )
    axes = axes.flatten()

    for ax, (fname, z) in zip(axes, subset):

        spec = load_spectrum(
            f"{base_path}/{fname}",
            z=z,
            **loader_kwargs
        )

        if output_flux_scale is None:
            output_flux_scale = spec.get("output_flux_scale")

        plot_spectrum_ax(
            ax,
            spec,
            title=fname,
            z=z,
            xlim=xlim,
            ylim=ylim,
            **plot_kwargs
        )

    for ax in axes[len(subset):]:
        ax.axis("off")

    fig.supxlabel(r"Rest-frame wavelength [$\mu$m]")
    if output_flux_scale is not None:
        power = -int(np.log10(output_flux_scale))
        ylabel = (
            rf"$F_\lambda$ "
            rf"($10^{{{power}}}$ erg s$^{{-1}}$ cm$^{{-2}}$ Å$^{{-1}}$)"
        )
    else:
        ylabel = r"$F_\lambda$ (erg s$^{-1}$ cm$^{-2}$ Å$^{-1}$)"

    fig.supylabel(ylabel, fontsize=14)

    #fig.supylabel(r"Flux [10$^{-20}$ erg s$^{-1}$ cm$^{-2}$ Å$^{-1}$]")

    plt.tight_layout()
    return fig
