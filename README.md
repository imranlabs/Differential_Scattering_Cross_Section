<!-- Open in Colab for your main notebook -->
  <a href="https://colab.research.google.com/github/imranlabs/Differential_Scattering_Cross_Section/blob/main/Differential_Scattering_Cross_Section.ipynb">
    <img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab">
  </a>
</p>

<p align="left">
  <!-- Python version -->
  <a href="https://www.python.org/">
    <img src="https://img.shields.io/badge/Python-3.8%2B-blue" alt="Python 3.8+">
  </a>
  <!-- License (auto-reads your repo’s license) -->
  <a href="https://github.com/imranlabs/Differential_Scattering_Cross_Section/blob/main/LICENSE">
    <img src="https://img.shields.io/github/license/imranlabs/Broadband_Plasmonic_cloaking" alt="License">
  </a>

This repository provides a **fully reproducible** Jupyter notebook and minimal helper functions to compute and visualize the **differential scattering cross section** (far‑field angular distribution) for a dielectric core surrounded by a **concentric, disordered shell** of gold nanoparticles (AuNPs). Using a multiscale computational framework, we quantify how controlled disorder in the AuNP shell redistributes scattered power, suppresses side lobes, and biases scattering forward. This tunability suggests potential applications in passive cloaking and high-resolution imaging. The workflow is geared toward quick parameter sweeps, clean plots, and exportable data for downstream analysis.

> **Paper**: Md. Imran Khan, Sayantani Ghosh, and Arnold D. Kim, *“Scattering by nanoplasmonic mesoscale assemblies,”* **JOSA A** 42, 1244–1253 (2025). DOI: https://doi.org/10.1364/JOSAA.560629

## Key Features

- Multiscale Modeling: Combines nanoparticle–light interactions (Foldy–Lax multiple scattering) and core response (Method of Fundamental Solutions).

Parameter Studies: Sweeps AuNP diameters (5–20 nm) and filling fractions (0.1–0.3) to compute angular scattering patterns.

- ### Physical Quantities Computed:

    - Differential scattering cross-section($\sigma_{d}$)

    - Total scattering cross-section($\sigma_{t}$)

    - Scattering albedo($\bar{\omega_{0}}$)

    - Anisotropy factor ($g$)

- **Henyey-Greenstein Fitting:** Used to extract anisotropy from scattering profiles.

## Repository Structure
```
├── Differential_scattering_cross_section.ipynb  # Main Jupyter notebook
├── src/
│   ├── Geometry.py                          # Geometry-related functions
│   └── Mfs.py                               # Method of Fundamental Solutions
├── data/
│   └── optical_constants/
│   |    └── 20nm_gold_film_silica.csv       # Refractive index data
|   └-- plot_data/
│       └── Core_750_5nm_20nm_in_670nm.csv   # far-field scattering data
|       └── Core_750_5nm_20nm_in_450nm.csv                 
|                    
├── results/
│   ├── data/                                # Output CSV files
│   └── figures/                             # output plots
└── README.md
```
## Dependencies

 - Python 3.x

 - NumPy

 - Matplotlib

 - Pandas

 - SciPy

 - tqdm

**Install missing packages using:**
`pip install numpy matplotlib pandas scipy tqdm`

## Usage

1. Open the Jupyter notebook:
    `Differential_scattering_cross_section.ipynb`
   

3. Run the cells sequentially:

 - Environment Setup: Imports necessary libraries.

 - Parameters/Config: Loads refractive index data and sets geometric/optical parameters.

 - Main Execution Loop: Computes core-only and core-shell scattering properties.

- Plotting: Generates wavelength-dependent anisotropy plots, superimposed polar plots, and on-demand polar plots.


3. Outputs:
   
The notebook highlights how **particle size**, **filling fraction/number density**, and **shell thickness** affect the **angular distribution** of scattered light. Typical outputs include:

 - CSV files with computed scattering properties (in results/data/)

 - Figures comparing anisotropy factors (in results/figures/)

## Citation
If you use this notebook or derived data/figures, please cite the article:

```
Md. Imran Khan, Sayantani Ghosh, and Arnold D. Kim,
"Scattering by nanoplasmonic mesoscale assemblies,"
J. Opt. Soc. Am. A 42, 1244–1253 (2025).
https://doi.org/10.1364/JOSAA.560629
```

## Acknowledgments

This work builds on earlier research:
- [Modeling broadband cloaking using 3D nano-assembled plasmonic meta-structures](https://github.com/imranlabs/Broadband_Plasmonic_cloaking)
- Many thanks to collaborators and the broader plasmonics community.

## License
This project is licensed under the MIT License. See LICENSE file for details.


