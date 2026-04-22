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

    # --- Fixed rest-frame emission lines (μm) ---
    emission_lines = {
        r"[O II]": 0.3727,
        r"[Ne III]": 0.386876,
        r"[O III] 4363": 0.436321,
        r"[O III] 5007": 0.5006843,   	
    }

    H_emission_lines = {
        r"H$\beta$": 0.48613,
        r"H$\alpha$": 0.6563,
        r"H$\epsilon$": 0.3970079,
        r"H$\gamma$": 0.4340471,
        r"H$\delta$": 0.4101742,
    }
    
    # --- Draw emission lines first (background) ---
    for label, wave0 in emission_lines.items():
        ax.axvline(
            wave0,
            color="black",
            ls="--",
            lw=0.8,
            alpha=0.7,
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

    # --- Draw H series emission lines first (background) ---
    for label, wave0 in H_emission_lines.items():
        ax.axvline(
            wave0,
            color="deeppink",
            ls="--",
            lw=0.8,
            alpha=0.7,
            zorder=0
        )

        ax.text(
            wave0,
            0.78,
            label,
            rotation=90,
            ha="right",
            va="top",
            transform=ax.get_xaxis_transform(),
            fontsize=8,
            color="deeppink"
        )

    ax.step(
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
        r"[Ne III]": 0.386876,
        r"[O III] 4363": 0.436321,
        r"[O III] 5007": 0.5006843,   	
    }

    H_emission_lines = {
        r"H$\beta$": 0.48613,
        r"H$\alpha$": 0.6563,
        r"H$\epsilon$": 0.3970079,
        r"H$\gamma$": 0.4340471,
        r"H$\delta$": 0.4101742,
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

    # --- Draw H series emission lines first (background) ---
    for label, wave0 in H_emission_lines.items():
        ax.axvline(
            wave0,
            color="deeppink",
            ls="--",
            lw=0.8,
            alpha=0.6,
            zorder=0
        )

        ax.text(
            wave0,
            0.80,
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
            y_offset = 8.5 * (i - indices[0])

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

def plot_spectrum_presentation(
    fname,
    z,
    base_path="DeGraaff_espectros",
    lines=None,
    loader_kwargs=None,
    figsize=(8,5),
    xlim=None,
    ylim=None,
    title=None,
    spectrum_kwargs=None,
    line_kwargs=None,
):
    """
    Plot a single spectrum optimized for presentations.

    Parameters
    ----------
    fname : str
        Spectrum filename
    z : float
        Redshift
    lines : dict or None
        Emission lines to display {label: wavelength}
        If None, a default set is used
    spectrum_kwargs : dict
        kwargs passed to ax.step()
    line_kwargs : dict
        kwargs passed to ax.axvline()
    """

    import os
    import matplotlib.pyplot as plt

    if loader_kwargs is None:
        loader_kwargs = {}

    if spectrum_kwargs is None:
        spectrum_kwargs = dict(color="black", lw=1.5, where="mid")

    if line_kwargs is None:
        line_kwargs = dict(color="gray", ls="--", lw=1, alpha=0.7)

    # Default emission lines
    if lines is None:
        lines = {
            r"Ly$\alpha$": 0.121567,
            r"N V": 0.1240,
            r"Si IV": 0.140277,
            r"C IV": 0.1549,
            r"He II": 0.1640,
            r"C III]": 0.1909,
            r"Mg II]": 0.2800,
            r"[O II]": 0.3727,
            r"[Ne III]": 0.386876,
            r"H$\epsilon$": 0.3970079,
            r"H$\delta$": 0.4101742,
            r"H$\gamma$": 0.4340471,
            r"H$\beta$": 0.48613,
            r"[O III]": 0.5006843,
            r"He I 5876": 0.5875624,
            r"[O I]": 0.6300,
            r"H$\alpha$": 0.6563,
            r"[S II]": 0.6723,
            r"He I 7065": 0.7065,
            r"O I 8446": 0.8446,
            r"He I 10030": 1.0030,
            r"Pa$\delta$": 1.0049,
            r"He I 10830": 1.0830,
            r"Pa$\gamma$": 1.0938,
        }

    full_path = os.path.join(base_path, fname)

    spec = load_spectrum(
        full_path,
        z=z,
        **loader_kwargs
    )

    fig, ax = plt.subplots(figsize=figsize)

    # Plot spectrum
    ax.step(
        spec["wave"],
        spec["flux"],
        **spectrum_kwargs
    )

    # ---- Sort lines by wavelength ----
    sorted_lines = sorted(lines.items(), key=lambda x: x[1])

    # Levels used to avoid overlap
    levels = [0.98, 0.78, 0.98, 0.78]
    ha_lines = ['right', 'left', 'right', 'left']

    prev_wave = None
    level_index = 0

    for label, wave0 in sorted_lines:

        # Detect if lines are very close
        if prev_wave is not None and abs(wave0 - prev_wave) < 0.017:
            level_index += 1
        else:
            level_index = 0

        y = levels[level_index % len(levels)]
        ha = ha_lines[level_index % len(ha_lines)]

        # --- choose color ---
        if "H$" in label:
            color = "red"
            text_color = "red"
        else:
            color = "gray"
            text_color = "black"

        # vertical line
        ax.axvline(
            wave0,
            color=color,
            ls="--",
            lw=1,
            alpha=0.6
        )   

        ax.text(
            wave0,
            y,
            label,
            rotation=90,
            ha=ha,
            va="top",
            transform=ax.get_xaxis_transform(),
            fontsize=11,
            color=text_color
        )

        prev_wave = wave0

    if xlim is not None:
        ax.set_xlim(xlim)

    if ylim is not None:
        ax.set_ylim(ylim)

    ax.set_xlabel(r"Rest-frame wavelength [$\mu$m]", fontsize=13)
    ax.set_ylabel(r"Normalized $F_\lambda$", fontsize=13)

    if title is None:
        title = short_label_from_filename(fname)

    ax.set_title(f"{title} (z = {z:.3f})", fontsize=14)

    ax.grid(alpha=0.3)

    # --- Tick configuration ---
    ax.minorticks_on()

    ax.tick_params(
        axis="both",
        which="major",
        direction="in",
        length=6,
        width=1.2,
        labelsize=11,
        top=True, right=True
    )

    ax.tick_params(
        axis="both",
        which="minor",
        direction="in",
        length=3,
        width=1,
        top=True, right=True
    )

    # --- Border thickness ---
    for spine in ax.spines.values():
        spine.set_linewidth(1.2)

    fig.tight_layout()

    return fig, ax

def plot_spectrum_shaded_lines(
    fname,
    z,
    base_path="DeGraaff_espectros",
    lines=None,
    loader_kwargs=None,
    figsize=(8,5),
    xlim=None,
    ylim=None,
    title=None,
    spectrum_kwargs=None,
    shade_width=0.002,
):
    """
    Plot spectrum with shaded regions marking emission lines.
    Hydrogen lines are shown in red, others in gray.
    """

    import os
    import matplotlib.pyplot as plt

    if loader_kwargs is None:
        loader_kwargs = {}

    if spectrum_kwargs is None:
        spectrum_kwargs = dict(color="black", lw=1.4, where="mid")

    # Default emission lines
    if lines is None:
        lines = {
            r"[O II]": 0.3727,
            r"[Ne III]": 0.386876,
            r"H$\epsilon$": 0.3970079,
            r"H$\delta$": 0.4101742,
            r"H$\gamma$": 0.4340471,
            r"H$\beta$": 0.48613,
            r"[O III]": 0.5006843,
            r"He I 5876": 0.5875624,
            r"[O I]": 0.6300,
            r"H$\alpha$": 0.6563,
            r"[S II]": 0.6723,
            r"He I 7065": 0.7065,
            r"He I 10830": 1.0830,
        }

    full_path = os.path.join(base_path, fname)

    spec = load_spectrum(
        full_path,
        z=z,
        **loader_kwargs
    )

    fig, ax = plt.subplots(figsize=figsize)

    # Plot spectrum
    ax.step(
        spec["wave"],
        spec["flux"],
        **spectrum_kwargs
    )

    # ---- Sort lines by wavelength ----
    sorted_lines = sorted(lines.items(), key=lambda x: x[1])

    # label vertical levels
    levels = [0.98, 0.88, 0.98, 0.88]

    prev_wave = None
    level_index = 0

    for label, wave0 in sorted_lines:

        if prev_wave is not None and abs(wave0 - prev_wave) < 0.015:
            level_index += 1
        else:
            level_index = 0

        y = levels[level_index % len(levels)]

        # choose color
        if "H$" in label:
            color = "firebrick"
        else:
            color = "gray"

        # shaded region
        ax.axvspan(
            wave0 - shade_width,
            wave0 + shade_width,
            color=color,
            alpha=0.18,
            zorder=0
        )

        # central line
        ax.axvline(
            wave0,
            color=color,
            lw=1,
            alpha=0.4,
            zorder=1
        )

        # label
        ax.text(
            wave0,
            y,
            label,
            rotation=90,
            ha="center",
            va="top",
            transform=ax.get_xaxis_transform(),
            fontsize=11,
            color=color
        )

        prev_wave = wave0

    if xlim is not None:
        ax.set_xlim(xlim)

    if ylim is not None:
        ax.set_ylim(ylim)

    # Axis labels
    ax.set_xlabel(r"Rest-frame wavelength [$\mu$m]", fontsize=13)
    ax.set_ylabel(r"Flux", fontsize=13)

    if title is None:
        title = short_label_from_filename(fname)

    ax.set_title(title, fontsize=14)

    ax.text(
        0.15, 0.95,
        f"z = {z:.3f}",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=12
    )

    # ---- Style adjustments ----
    ax.minorticks_on()

    ax.tick_params(
        axis="both",
        which="major",
        direction="in",
        length=6,
        width=1.2,
        labelsize=11,
        top=True,
        right=True
    )

    ax.tick_params(
        axis="both",
        which="minor",
        direction="in",
        length=3,
        width=1
    )

    # thicker border
    for spine in ax.spines.values():
        spine.set_linewidth(1.2)

    ax.grid(alpha=0.25)

    fig.tight_layout()

    return fig, ax

def plot_mean_spectrum(
    mean_spec,
    lines=None,
    figsize=(8,5),
    xlim=None,
    ylim=None,
    title=None,
    spectrum_kwargs=None,
    std_kwargs=None,
    err_kwargs=None,
    line_kwargs=None,
    min_contrib=None,
    show_std=True,    
    show_err=False,
):
    """
    Plot mean spectrum with optional std shading.

    Parameters
    ----------
    ax : matplotlib axis
    mean_spec : dict
        Output of compute_mean_spectrum()
    show_std : bool
        Plot shaded std region
    show_err : bool
        Plot error on the mean
    min_contrib : int or None
        Mask regions with fewer contributing spectra
    lines : dict or None
        Emission lines {label: wavelength}
        None → usa default
        {} → não plota linhas
    """

    wave = mean_spec["wave"]
    flux = mean_spec["flux_mean"]
    std  = mean_spec.get("flux_std")
    err  = mean_spec.get("err_mean")
    n_contrib = mean_spec.get("n_contrib")

    # -------------------------
    # kwargs defaults
    # -------------------------
    if spectrum_kwargs is None:
        spectrum_kwargs = dict(color="black", lw=1.5, where="mid")

    if std_kwargs is None:
        std_kwargs = dict(color="gray", alpha=0.3)

    if err_kwargs is None:
        err_kwargs = dict(color="red", alpha=0.2)

    if line_kwargs is None:
        line_kwargs = dict(ls="--", lw=1, alpha=0.6)

    # -------------------------
    # máscara de qualidade
    # -------------------------
    if min_contrib is not None and n_contrib is not None:
        mask = n_contrib >= min_contrib
        wave = wave[mask]
        flux = flux[mask]
        if std is not None:
            std = std[mask]
        if err is not None:
            err = err[mask]

    # -------------------------
    # default lines
    # -------------------------
    if lines is None:
        lines = {
            r"[O II]": 0.3727,
            r"[Ne III]": 0.386876,
            r"H$\epsilon$": 0.3970079,
            r"H$\delta$": 0.4101742,
            r"H$\gamma$": 0.4340471,
            r"H$\beta$": 0.48613,
            r"[O III]": 0.5006843,
            r"H$\alpha$": 0.6563,
        }

    fig, ax = plt.subplots(figsize=figsize)

    # -------------------------
    # plot espectro médio
    # -------------------------
    ax.step(wave, flux, **spectrum_kwargs)

    # -------------------------
    # std shading
    # -------------------------
    if show_std and std is not None:
        ax.fill_between(
            wave,
            flux - std,
            flux + std,
            **std_kwargs
        )

    # -------------------------
    # erro da média
    # -------------------------
    if show_err and err is not None:
        ax.fill_between(
            wave,
            flux - err,
            flux + err,
            **err_kwargs
        )

    # -------------------------
    # linhas de emissão
    # -------------------------
    if lines:  # vazio {} → não entra
        sorted_lines = sorted(lines.items(), key=lambda x: x[1])

        levels = [0.98, 0.78, 0.98, 0.78]
        ha_lines = ['right', 'left', 'right', 'left']

        prev_wave = None
        level_index = 0

        for label, wave0 in sorted_lines:

            if prev_wave is not None and abs(wave0 - prev_wave) < 0.017:
                level_index += 1
            else:
                level_index = 0

            y = levels[level_index % len(levels)]
            ha = ha_lines[level_index % len(ha_lines)]

            # cores
            if "H$" in label:
                color = "red"
                text_color = "red"
            else:
                color = "gray"
                text_color = "black"

            ax.axvline(wave0, color=color, **line_kwargs)

            ax.text(
                wave0,
                y,
                label,
                rotation=90,
                ha=ha,
                va="top",
                transform=ax.get_xaxis_transform(),
                fontsize=11,
                color=text_color
            )

            prev_wave = wave0

    # -------------------------
    # estética
    # -------------------------
    if xlim is not None:
        ax.set_xlim(xlim)

    if ylim is not None:
        ax.set_ylim(ylim)

    ax.set_xlabel(r"Rest-frame wavelength [$\mu$m]", fontsize=13)
    ax.set_ylabel(r"Mean normalized $F_\lambda$", fontsize=13)

    if title is not None:
        ax.set_title(title, fontsize=14)

    ax.grid(alpha=0.3)

    ax.minorticks_on()

    ax.tick_params(
        axis="both",
        which="major",
        direction="in",
        length=6,
        width=1.2,
        labelsize=11,
        top=True, right=True
    )

    ax.tick_params(
        axis="both",
        which="minor",
        direction="in",
        length=3,
        width=1,
        top=True, right=True
    )

    for spine in ax.spines.values():
        spine.set_linewidth(1.2)

    if "n_objects" in mean_spec:
        ax.text(
            0.85, 0.95,
            f"N = {mean_spec['n_objects']}",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=11,
            bbox=dict(facecolor='white', alpha=0.7, edgecolor='none')
        )

    fig.tight_layout()

    return fig, ax

def plot_overlaid_mean_spectra(
    mean_specs,
    xlim=(0.2, 0.6),
    ylim=None,
    figsize=(7, 5),
    offset=True,
    lines=None,
    cmap_name="RdPu_r",
    min_contrib=None,
):
    """
    Plot overlaid mean spectra for multiple groups.

    Parameters
    ----------
    mean_specs : dict
        {"G1": mean_spec, "G2": mean_spec, ...}
    offset : bool
        Apply vertical offset to separate spectra
    lines : dict or None
        Emission lines {label: wavelength}
    """

    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib import cm

    # -------------------------
    # cores
    # -------------------------
    n = len(mean_specs)
    cmap = cm.get_cmap(cmap_name)
    colors = cmap(np.linspace(0.0, 0.8, n))

    fig, ax = plt.subplots(figsize=figsize)

    ax.set_xlim(*xlim)
    if ylim is not None:
        ax.set_ylim(*ylim)

    xmin, xmax = ax.get_xlim()

    # -------------------------
    # linhas de emissão
    # -------------------------
    if lines is None:
        lines = {
            r"[O II]": 0.3727,
            r"[Ne III]": 0.386876,
            r"H$\delta$": 0.4101742,
            r"H$\gamma$": 0.4340471,
            r"H$\beta$": 0.48613,
            r"[O III]": 0.5006843,
        }

    if lines:
        for label, wave0 in lines.items():
            ax.axvline(wave0, color="gray", ls="--", lw=0.8, alpha=0.6, zorder=0)

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

    # -------------------------
    # plot dos grupos
    # -------------------------
    for j, (name, mean_spec) in enumerate(mean_specs.items()):

        wave = mean_spec["wave"]
        flux = mean_spec["flux_mean"]
        n_contrib = mean_spec.get("n_contrib")

        # máscara de qualidade
        if min_contrib is not None and n_contrib is not None:
            mask = n_contrib >= min_contrib
            wave = wave[mask]
            flux = flux[mask]

        color = colors[j]

        if offset:
            y_offset = j * 1.5  # ajuste fino aqui
            y = flux + y_offset

            ax.step(wave, y, where="mid", color=color, lw=1.5)

            # posição do label
            axis_width = xmax - xmin
            free_space = xmax - wave.max()

            if free_space >= 0.2 * axis_width:
                x_text = wave.max() * 1.01
                ha = "left"
            else:
                x_text = wave.max() * 0.99
                ha = "right"

            y_text = np.nanmedian(y[-50:])  # posição robusta

            ax.text(
                x_text,
                y_text,
                f"{name} (N={mean_spec['n_objects']})",
                fontsize=9,
                ha=ha,
                va="center"
            )

        else:
            ax.step(
                wave,
                flux,
                where="mid",
                color=color,
                lw=1.8,
                label=f"{name} (N={mean_spec['n_objects']})"
            )

    # -------------------------
    # estética
    # -------------------------
    ax.set_xlim(*xlim)

    if ylim is not None:
        ax.set_ylim(*ylim)

    ax.set_xlabel(r"Rest-frame wavelength [$\mu$m]")
    ax.set_ylabel(r"Normalized Flux")

    if not offset:
        ax.legend(frameon=False)

    fig.tight_layout()

    return fig, ax