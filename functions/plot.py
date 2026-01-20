import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from .spectrum import load_spectrum
import os


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
    Plot a panel of spectra.

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

    n_panel = nrows * ncols
    subset = spec_info[start:start + n_panel]

    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(14, 4 * nrows),
        sharex=True, sharey=True
    )
    axes = axes.flatten()

    output_flux_scale = None
    normalized = None
    norm_window = None

    for ax, (fname, z) in zip(axes, subset):

        spec = load_spectrum(
            f"{base_path}/{fname}",
            z=z,
            **loader_kwargs
        )

        # ---- Skip spectra that failed normalization ----
        if loader_kwargs.get("normalize", False) and not spec["normalized"]:
            print(
                f"⚠️ Skipping {fname}, redshift {z}: normalization failed "
                f"({spec['norm_error']})"
            )
            ax.axis("off")
            continue

        # ---- Consistency checks ----
        if normalized is None:
            normalized = spec.get("normalized", False)
            output_flux_scale = spec.get("output_flux_scale")
            norm_window = spec.get("norm_window")
        else:
            if spec.get("normalized", False) != normalized:
                raise ValueError("Inconsistent normalization across spectra")
            if spec.get("output_flux_scale") != output_flux_scale:
                raise ValueError("Inconsistent output_flux_scale across spectra")

        plot_spectrum_ax(
            ax,
            spec,
            title=fname,
            z=z,
            xlim=xlim,
            ylim=ylim,
            **plot_kwargs
        )
    
    # Turn off unused axes
    for ax in axes[len(subset):]:
        ax.axis("off")

    # ---- Axis labels ----
    fig.supxlabel(r"Rest-frame wavelength [$\mu$m]")

    if normalized:
        ylabel = r"Normalized $F_\lambda$"
    elif output_flux_scale is not None:
        power = -int(np.log10(output_flux_scale))
        ylabel = (
            rf"$F_\lambda$ "
            rf"($10^{{{power}}}$ erg s$^{{-1}}$ cm$^{{-2}}$ Å$^{{-1}}$)"
        )
    else:
        ylabel = r"$F_\lambda$ (erg s$^{-1}$ cm$^{-2}$ Å$^{-1}$)"

    fig.supylabel(ylabel, fontsize=14)


    plt.tight_layout()
    return fig

def short_label_from_filename(fname):
    """
    Convert:
    capers-cos07-v4_prism-clear_6368_43711.spec.fits
    → capers-6368_43711
    """

    base = os.path.basename(fname)

    # remove sufixo
    base = base.replace(".spec.fits", "")

    # survey = primeira parte antes do primeiro "-"
    survey = base.split("-")[0]

    # pega os dois últimos blocos separados por "_"
    id_part = "_".join(base.split("_")[-2:])

    return f"{survey}-{id_part}"


def plot_overlaid_spectra(
    spec_info,
    indices,
    base_path="DeGraaff_espectros",
    loader_kwargs=None,
    xlim=(0.2, 0.6),
    ylim=None,
    figsize=(7, 5),
    offset=True,
):
    if loader_kwargs is None:
        loader_kwargs = {}

     # --- Fixed rest-frame emission lines (μm) ---
    emission_lines = {
        r"[O II]": 0.3727,
        r"[O III] 4363": 0.4363,
        r"[O III] 5007": 0.5007,
        r"H$\beta$": 0.4861,
        r"H$\alpha$": 0.6563,
    }
    
    n = len(indices)
    base_cmap = cm.get_cmap("RdPu_r")  # escolha a paleta aqui
    colors = base_cmap(np.linspace(0.0, 0.7, n))


    fig, ax = plt.subplots(figsize=figsize)

    # --- FIXO: definir limites ANTES de usar get_xlim/get_ylim ---
    ax.set_xlim(*xlim)
    if ylim is not None:
        ax.set_ylim(*ylim)

    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()

    # --- Draw emission lines first (background) ---
    for label, wave0 in emission_lines.items():
        ax.axvline(
            wave0,
            color="gray",
            ls="--",
            lw=0.8,
            alpha=0.6,
            zorder=0
        )

        ax.text(
            wave0,
            0.98,
            label,
            rotation=90,
            ha="right",
            va="top",
            transform=ax.get_xaxis_transform(),
            fontsize=8,
            color="gray"
        )

    # --- Plot spectra ---
    for j, i in enumerate(indices):
        fname, z = spec_info[i]

        color = colors[j]  # cor específica para este espectro

        # Convert numpy string → Python string
        fname = str(fname)

        full_path = os.path.join(base_path, fname)

        try:
            data = load_spectrum(
                full_path,
                z=z,
                **loader_kwargs
            )
        except Exception as e:
            print(f"Skipping {fname} → {e}")
            continue

        if offset:
            label = short_label_from_filename(fname)
            y_offset = 2.5 * (i - indices[0])

            y = data["flux"] + y_offset
            x = data["wave"]

            #ax.plot(x, y, lw=1)
            ax.step(x, y, color=color, where="mid", lw=1)

            # posição do texto (lado direito do gráfico)
            #xmin, xmax = ax.get_xlim()

            axis_width = xmax - xmin
            free_space = xmax - x.max()

            if free_space >= 0.2 * axis_width: 
                # cabe texto fora do espectro
                x_text = x.max() * 1.01
                ha = "left"
            else:
                # coloca dentro do espectro
                x_text = x.max() * 0.99
                ha = "right"

            dy = 1.15
            y_text = np.nanmedian(y[-100:-5]) + dy

            ax.text(
                x_text,
                y_text,
                label,
                fontsize=8,
                ha=ha,
                va="center"
            )

        else:
            label = short_label_from_filename(fname)
            ax.step(data["wave"],data["flux"],where='mid',lw=1.2,label=label)

    ax.set_xlim(*xlim)
    if ylim is not None:
        ax.set_ylim(*ylim)

    ax.set_xlabel(r"Rest-frame wavelength [$\mu$m]")
    ax.set_ylabel(r"Normalized Flux (arbitrary units)")

    ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()

    return fig
