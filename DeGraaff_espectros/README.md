This repository is a snapshot of the public JWST NIRSpec spectra processed with the msaexp pipeline. Please refer to and cite de Graaff et al. (2024) and Heintz et al. (2025) for the main presentation of the msaexp pipeline.

This release corresponds to the "v4" version of the spectral extractions, which significantly extends the wavelength range of the extracted spectra to regions that may suffer contamination of overlapping spectral orders.  The sensitivity of the higher orders is strongly weighted toward the blue side of the spectrum for all gratings, so, in practice, relatively red galaxies often suffer relatively minor order contamination.

Please refer to and cite this [DOI 10.5281/zenodo.1547235](https://zenodo.org/records/15472354) and [Valentino et al. (2025)](https://ui.adsabs.harvard.edu/abs/2025A&A...699A.358V) when using this specific data release and the ``v4`` spectra for a presentation of the extended extractions. For more information and updates, please refer to the DJA Blog Post: https://dawn-cph.github.io/dja/blog/2025/05/01/nirspec-merged-table-v4/



## Data Content

This release provides version 4 of the NIRSpec Merged Table, a comprehensive catalog of uniformly reduced and analyzed JWST/NIRSpec spectra.
The data have been processed using the [msaexp](https://github.com/gbrammer/msaexp) and [grizli](https://github.com/gbrammer/grizli) pipelines and are publicly available through the DAWN JWST Archive (DJA).
The catalog integrates multiple data products, including spectral extractions, redshift measurements, emission line fluxes, and photometric associations.

The merged table consolidates information from several database tables:

- `nirspec_extractions`: Basic spectrum parameters (e.g., grating, mask, exposure time).
- `nirspec_redshifts`: Redshift fit results and emission line fluxes.
- `nirspec_redshifts_manual`: Grades and comments from visual inspection.
- `nirspec_integrated`: Observed- and rest-frame filters integrated through the spectra at the derived redshift.
- `grizli_photometry`: Photometry and some EAZY outputs of the nearest counterpart in the DJA/grizli photometric catalogs.

The catalog includes approximately 80,367 entries, each corresponding to a unique NIRSpec spectrum.
Each entry contains metadata such as source ID, coordinates, grating/filter configuration, exposure time, redshift estimates, emission line measurements, and photometric associations.

```
| N     |  Grating-Filter | Concatenated 1D spectrum file |
|------:|:---------------:|:------------------------------|
|   113 |  G140H-F070LP   |  dja_msaexp_emission_lines_v4.4.g140h-f070lp_spectra.fits  |
|   684 |  G140H-F100LP   |  dja_msaexp_emission_lines_v4.4.g140h-f100lp_spectra.fits  |
|  5851 |  G140M-F070LP   |  dja_msaexp_emission_lines_v4.4.g140m-f070lp_spectra.fits  |
|  2165 |  G140M-F100LP   |  dja_msaexp_emission_lines_v4.4.g140m-f100lp_spectra.fits  |
|  6179 |  G235H-F170LP   |  dja_msaexp_emission_lines_v4.4.g235h-f170lp_spectra.fits  |
|  8000 |  G235M-F170LP   |  dja_msaexp_emission_lines_v4.4.g235m-f170lp_spectra.fits  |
|  8820 |  G395H-F290LP   |  dja_msaexp_emission_lines_v4.4.g395h-f290lp_spectra.fits  |
| 13606 |  G395M-F290LP   |  dja_msaexp_emission_lines_v4.4.g395m-f290lp_spectra.fits  |
| 34949 |  PRISM-CLEAR    |  dja_msaexp_emission_lines_v4.4.prism_spectra.fits         |
```

A searchable HTML overview of the public spectra is available at [https://s3.amazonaws.com/msaexp-nirspec/extractions/nirspec_public_v4.4.html](https://s3.amazonaws.com/msaexp-nirspec/extractions/nirspec_public_v4.4.html).


## Usage Notes

- Data Format: The main catalog is provided in compressed CSV format. Column descriptions, including units and formats, are detailed in the accompanying columns CSV file.

The redshift and line fluxes are derived with the functions msaexp.fit_redshift and msaexp.plot_spectrum. The continuum is modelled as a combination of splines, whose coefficient are listed in the table. Emission lines are superimposed as Gaussian profiles, smoothed by the instrumental resolution (increased by a factor of 1.3x compared with the nominal line spread functions in the JWST User Documentation, JDox, de Graaff+2025) and a fixed line velocity width of 100 km/s. Equivalent widths in angstrom in the observed frame are also reported. A dictionary with the available lines and their rest-frame wavelengths in vacuum can be generated with grizli.utils.get_line_wavelengths(). The fit is performed with a least square template method. The uncertainties are rescaled with a polynomial curve prior to fitting such that '(flux - model)/(err*scl)' residuals are 'N(0,1)' using msaexp.spectrum.calc_uncertainty_scale. A systematic uncertainty floor is introduced. The robustness of the redshift solution is flagged according to the following scheme:
- **Grade 3**: Robust redshift from one or more emission absorption features
- **Grade 2**: Ambiguous continuum features, perhaps only one line or low confidence lines
- **Grade 1**: No clear features in the spectrum to constrain the redshift
- **Grade 0**: Spectrum suffers some data quality issue and should
- **Grade -1**: Fit not performed or graded

If multiple spectra of the same sources are available, a common best redshift solution is stored in the z_best column.

- Spectral Data: Each entry includes links to the corresponding 2D spectra in fits and png format, which can be accessed through the DJA interface or the public spectra overview page at: https://s3.amazonaws.com/msaexp-nirspec/extractions/nirspec_public_v4.4.html
The 1D spectra in the fits tables in this release are in the format generated by msaexp. The columns are described in msaexp_spectrum.columns.csv.

- If existing, a comparison with the previous v3 version the spectra generated with msaexp is available at the same link on the DJA interface.  

- Photometric Associations: Photometric data are matched to the nearest counterparts in the DJA/grizli catalogs, providing additional context for each spectroscopic observation. The latest public versions of the DJA/grizli catalogs can are described here: https://dawn-cph.github.io/dja/imaging/v7/
The latest v7/ version of the imaging mosaics and catalogs can be retrieved here: https://s3.amazonaws.com/grizli-v2/JwstMosaics/v7/index.html
A description of the data reduction and catalog creation is available on the DJA interface and in Valentino+2023.

## Caveats

- An effective extended-source path-loss correction for light outside of the slitlet for each source using the a priori position within the shutter and assuming an azimuthally symmetric Gaussian profile is applied to each spectrum in this release (de Graaff et al. 2025, Section 3). However, the spectra are not rescaled to match the observed photometry in the "phot_" columns.
- The flux calibration is derived from calibration, monitoring, and scientific programs (Valentino et al. 2025). However, residual features in spectra due to an imperfect cross-calibration of different overlapping orders are still present in the released spectra.
- Examples of second order corrections are presented in Ito et al. (2025), where medium-resolution grating spectra beyond their nominal coverage deviate from the prism counterpart (Figure C.1 in appendix in Ito et al. 2025). Moreover, a significant downturn is present in prism spectra at wavelengths longer than 5.2$\mu$. These second order corrections are largely mitigated cross-calibrating all the available spectra and photometry by means of simple low-order polynomial corrections, for example available in spectrophotometric modeling codes. Future calibration programs dedicated to reconstruction of the sensitivity curves in the extended spectra at different location of the MSA will allow for refined absolute calibrations.

## Demos

The steps described in the Usage Notes are collected in publicly available Jupyter notebook on DJA's website:

- DJA spectroscopic data products: https://dawn-cph.github.io/dja/blog/2023/07/18/nirspec-data-products/
- v4 Merged Table properties: https://dawn-cph.github.io/dja/blog/2025/05/01/nirspec-merged-table-v4/



## Software and Dependencies

The data processing and analysis utilized the following software packages:

- msaexp: v0.9.8.dev3+ge0e3f39.d20250429
- grizli: v1.12.12.dev5+g5896d62.d20250426
- eazy-py: v0.8.5


## Contact

For questions or feedback regarding this data release, please contact Gabriel Brammer at gabriel.brammer@nbi.ku.dk.


## Citations

We warmly encourage the users to refer to the original works describing the datasets collected in this release. Relevant references, compiled to the best of our knowledge, are available here: [LINK].

Works that make use of the products of the DJA should cite this DOI and relevant articles describing them:

1) Msaexp pipeline and methods, v3 and previous spectroscopic compilation releases:
- de Graaff, A., Brammer, G., Weibel, A.,  et al., "RUBIES: a complete census of the bright and red distant Universe with JWST/NIRSpec", A&A, 697, 189 (2025)
- Heintz, K. E., Brammer, G., Watson, D., et al., "The JWST-PRIMAL archival survey: A JWST/NIRSpec reference sample for the physical properties and Lyman-α absorption and emission of ∼600 galaxies at z = 5.0−13.4", A&A, 693, 60 (2025)
- Brammer G., "msaexp: NIRSpec analyis tools", 10.5281/zenodo.7299500 (2022)


2) Spectroscopic release v4 and extended spectra:
- Valentino, F., Heintz, K. E., Brammer, G. et al., "Gas outflows in two recently quenched galaxies at z = 4 and 7", A&A, 699, 358 (2025)
- Pollock, C., Gottumukkala, R., Heintz, K. E. et al., "Novel z~10 auroral line measurements extend the gradual offset of the FMR deep into the first Gyr of cosmic time ", arXiv:2506.15779 (2025)

3) Grizli pipeline:
- Brammer G., "grizli", https://zenodo.org/records/8370018 (2023)

4) Imaging release:
- Valentino, F., Brammer, G., Gould, K. M. L. et al., "An Atlas of Color-selected Quiescent Galaxies at z > 3 in Public JWST Fields", ApJ, 947, 20 (2023) 
    