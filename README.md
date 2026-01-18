# Columbus Rotational Scan Data — CHIME/FRB Catalog 2

## Overview

This repository contains the full reproducible workflow and data products for the paper:

**“Independent Confirmation of the Columbus Axis Using CHIME/FRB Catalog 2”**  
Zenodo DOI: https://doi.org/10.5281/zenodo.18286776

This work represents the third stage in a research program exploring weak directional structure in cosmological observables around a theoretically motivated CPT-symmetric axis — the **Siamese Axis**. After two previous detections using FRB–QSO correlations (v1.3) and CHIME/FRB Catalog 1 (v1.7), this repository provides an **independent confirmation** using **CHIME/FRB Catalog 2**, a larger and more heterogeneous dataset.

Using the **Columbus Rotational Scan**, a rotational hemispheric test developed in earlier studies, we detect a clean and reproducible sinusoidal modulation in FRB dispersion measures ΔDM(ψ), with:

- **R² = 0.87**  
- **Effective amplitude A ≈ 56.7 pc cm⁻³**  
- **Phase φ₀ ≈ 133.6°**, matching Catalog 1 and FRB–QSO results  
- **Permutation test output included**

The stability of the *phase* across three independent analyses strongly suggests the presence of a real, narrow-band directional structure in the sky.

---

## Repository Structure

## How to Reproduce the Analysis

1. **Create the project folder structure**, for example:
chime2/
data/
src/
figs/

*(The catalog cannot be redistributed and must be downloaded manually.)*

Place the analysis script in:chime2/src/analysis.py

Open **PowerShell** inside the `chime2` directory and install the required Python packages:

```powershell
pip install numpy pandas matplotlib scipy

The script will automatically:

load the catalog from data/

apply extragalactic filters (|b| > 20°, DM > 800)

compute the rotational angle φ around the Siamese Axis

perform the Columbus Rotational Scan (Mode B)

generate figures into figs/

generate numerical results into data/

After execution, you will find:

chime2/data/chime2_columbus_scan.csv

chime2/data/chime2_sinefit_summary.json

chime2/data/chime2_perm_amplitudes.npy

all figures in chime2/figs/


## Scientific Context and Interpretation

The Columbus Rotational Scan performed on CHIME/FRB Catalog 2 builds upon two previous studies:  
(1) *Testing CPT-Symmetric Siamese Universes through FRB–QSO Sky Correlations* (v1.3), and  
(2) *Rotational Hemispheric Test around a Siamese CPT-Symmetric Axis* (v1.7).  
Together, these works form a consistent sequence of independent analyses pointing toward the same preferred direction in the sky.

In the Siamese Cosmology framework, this directional structure is interpreted as the observational imprint of a phase-desynchronization mechanism Δφ(a) between two CPT-reflected cosmic sectors. The model predicts that such a mechanism should manifest not as a simple dipole, but as a **rotational (azimuthal) modulation** around a fixed axis at approximately (RA ≈ 170°, Dec ≈ 40°).

The results from CHIME/FRB Catalog 2 — specifically the sinusoidal ΔDM(ψ) modulation with phase φ₀ ≈ 133.6° and a high-quality fit (R² ≈ 0.87) — are **consistent with this theoretical prediction**. Importantly, the recovered phase aligns closely with the values obtained from both the FRB–QSO correlation study and CHIME/FRB Catalog 1.

While the Columbus Scan confirms the **phenomenon** (a stable, reproducible azimuthal anisotropy), the underlying **cause** — phase desynchronization — remains an interpretative element of the theoretical framework. Future datasets (CHIME/FRB Catalog 3, DSA-2000, SKA) will be required to test the mechanism more directly.

In summary, the Siamese Cosmology framework gains **empirical support** from the directional consistency observed across three independent analyses, but the physical origin of the anisotropy remains an open and testable question.

