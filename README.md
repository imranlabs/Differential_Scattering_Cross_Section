<p align="left">
  <!-- Python version -->
  <a href="https://www.python.org/">
    <img src="https://img.shields.io/badge/Python-3.8%2B-blue" alt="Python 3.8+">
  </a>
  <!-- License (auto-reads your repo’s license) -->
  <a href="https://github.com/imranlabs/Differential_Scattering_Cross_Section/blob/main/LICENSE">
    <img src="https://img.shields.io/github/license/imranlabs/Broadband_Plasmonic_cloaking" alt="License">
  </a>
  <!-- Open in Colab for your main notebook -->
  <a href="https://colab.research.google.com/github/imranlabs/Differential_Scattering_Cross_Section/blob/main/Differential_Scattering_Cross_Section.ipynb">
    <img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab">
  </a>
</p>

# Scattering by nanoplasmonic mesoscale assemblies: [Tuining Diffraction](https://doi.org/10.1364/JOSAA.560629)
Reproducible notebook for computing and visualizing the **differential scattering cross section(diffraction)** of a dielectric core coated with a disordered shell of gold nanoparticles (AuNPs).
This repository contains the code(`Differential_Scattering_Cross_Section.ipynb`) and data for the study "Scattering by nanoplasmonic mesoscale assemblies" [published in JOSA A (2025).](https://doi.org/10.1364/JOSAA.560629)

## Overview
We investigate the scattering behavior of nanoplasmonic mesoscale structures composed of a dielectric spherical core coated with a concentric shell containing randomly distributed gold nanoparticles (AuNPs). Using a multiscale computational framework, we quantify how controlled disorder in the AuNP shell redistributes scattered power, suppresses side lobes, and biases scattering forward. This tunability suggests potential applications in passive cloaking and high-resolution imaging.

## Key Features

- Multiscale Modeling: Combines nanoparticle–light interactions (Foldy–Lax multiple scattering) and core response (Method of Fundamental Solutions).

Parameter Studies: Sweeps AuNP diameters (5–20 nm) and filling fractions (0.1–0.3) to compute angular scattering patterns.

- ### Physical Quantities Computed:

    - Differential scattering cross-section($\sigma_{d}$)

    - Total scattering cross-section($\sigma_{t}$)

    - Scattering albedo($\bar{\omega_{0}}$)

    - Anisotropy factor (g)

- **Henyey-Greenstein Fitting:** Used to extract anisotropy from scattering profiles.

## Repository Structure
```none
├── Differential_scattering_cross_section.ipynb  # Main Jupyter notebook
├── src/
│   ├── Geometry.py                              # Geometry-related functions
│   └── Mfs.py                                  # Method of Fundamental Solutions
├── data/
│   └── optical_constants/
│   |    └── 20nm_gold_film_silica.csv          # Refractive index data
|   └-- plot_data/
│       └── Core_750_5nm_20nm_in_670nm.csv      # Pre-computed far-filed scattering data
|       └── Core_750_5nm_20nm_in_450nm.csv      # Pre-computed far-filed scattering data            
|                    
├── results/
│   ├── data/                                   # Output CSV files
│   └── figures/                                # Generated plots
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
    jupyter notebook Differential_scattering_cross_section.ipynb

2. Run the cells sequentially:

 - Environment Setup: Imports necessary libraries.

 - Parameters/Config: Loads refractive index data and sets geometric/optical parameters.

 - Main Execution Loop: Computes core-only and core-shell scattering properties.

- Plotting: Generates wavelength-dependent anisotropy plots, superimposed polar plots, on-demand polar plot.

3. Outputs:
 - CSV files with computed scattering properties (in results/data/)

 - Figures comparing anisotropy factors (in results/figures/)

## Citation

if you use this code or data, please cite:

 Imran Khan et al., “Mesoscale optical scattering control using nano-assembled plasmonic shells,” JOSA A (2025).
 DOI: https://doi.org/10.1364/JOSAA.560629

## Acknowledgments

This work builds on earlier research:
- [Modeling broadband cloaking using 3D nano-assembled plasmonic meta-structures](https://doi.org/10.1364/OE.395840)
- [Simulation](https://github.com/imranlabs/Broadband_Plasmonic_cloaking)

## License
This project is licensed under the MIT License. See LICENSE file for details.


