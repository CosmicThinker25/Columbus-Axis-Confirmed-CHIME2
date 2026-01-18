
# Columbus Rotational Scan Data — CHIME/FRB Catalog 2

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18286776.svg)](https://doi.org/10.5281/zenodo.18286776)
[![Science](https://img.shields.io/badge/Science-Open_Access-green.svg)](https://github.com/CosmicThinker25/Columbus-Axis-Confirmed-CHIME2)
[![Python](https://img.shields.io/badge/Python-3.9%2B-blue.svg)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## 🌌 Overview

This repository contains the full reproducible workflow and data products for the paper:

> **“Independent Confirmation of the Columbus Axis Using CHIME/FRB Catalog 2”**

This work represents the third stage in a research program exploring weak directional structure in cosmological observables around a theoretically motivated CPT-symmetric axis — the **Siamese Axis**. After two previous detections using FRB–QSO correlations (v1.3) and CHIME/FRB Catalog 1 (v1.7), this repository provides an **independent confirmation** using **CHIME/FRB Catalog 2**, a larger and more heterogeneous dataset.

Using the **Columbus Rotational Scan**, a rotational hemispheric test developed in earlier studies, we detect a clean and reproducible sinusoidal modulation in FRB dispersion measures $\Delta DM(\psi)$, with:

* **$R^2 = 0.87$** (High quality fit)
* **Effective amplitude $A \approx 56.7$ pc cm⁻³**
* **Phase $\phi_0 \approx 133.6^\circ$**, matching Catalog 1 and FRB–QSO results
* **Permutation test output included**

The stability of the *phase* across three independent analyses strongly suggests the presence of a real, narrow-band directional structure in the sky.

---

## 📂 Repository Structure

The project is organized as follows:

```text
.
├── data/
│   ├── chime2_columbus_scan.csv       # Output: Scan results
│   ├── chime2_sinefit_summary.json    # Output: Fit parameters
│   └── chime2_perm_amplitudes.npy     # Output: Permutation test data
├── src/
│   └── analysis.py                    # Main analysis script
├── figs/
│   └── (Generated figures will appear here)

Note: The original CHIME/FRB Catalog 2 file (chimefrbcat2.csv) cannot be redistributed here due to size/license and must be downloaded manually (see instructions below).

How to Reproduce the Analysis

Follow these steps to replicate the findings on your local machine.
1. Environment Setup

Open your terminal (PowerShell or Bash) and install the required scientific libraries:

pip install numpy pandas matplotlib scipy

2. Prepare the Data

Since the catalog file is external, you must:

    Download the CHIME/FRB Catalog 2 (CSV format) from the official website.

    Place the file chimefrbcat2.csv inside the data/ folder.

3. Run the Columbus Scan

Execute the main script from the project root:

python src/analysis.py

4. What Happens Next?

The script will automatically:

    Load the catalog from data/.

    Apply extragalactic filters (∣b∣>20∘, DM>800).

    Compute the rotational angle ψ around the Siamese Axis.

    Perform the Columbus Rotational Scan (Mode B).

    Generate figures in figs/ and numerical results in data/.

Scientific Context and Interpretation

The Columbus Rotational Scan performed on CHIME/FRB Catalog 2 builds upon two previous studies:

    Testing CPT-Symmetric Siamese Universes through FRB–QSO Sky Correlations (v1.3)

    Rotational Hemispheric Test around a Siamese CPT-Symmetric Axis (v1.7)

Together, these works form a consistent sequence of independent analyses pointing toward the same preferred direction in the sky.
The Theoretical Framework

In the Siamese Cosmology framework, this directional structure is interpreted as the observational imprint of a phase-desynchronization mechanism Δϕ(a) between two CPT-reflected cosmic sectors. The model predicts that such a mechanism should manifest not as a simple dipole, but as a rotational (azimuthal) modulation around a fixed axis at approximately (RA≈170∘, Dec≈40∘).
Consistency of Results

The results from CHIME/FRB Catalog 2 — specifically the sinusoidal ΔDM(ψ) modulation with phase ϕ0​≈133.6∘ and a high-quality fit (R2≈0.87) — are consistent with this theoretical prediction. Importantly, the recovered phase aligns closely with the values obtained from both the FRB–QSO correlation study and CHIME/FRB Catalog 1.

While the Columbus Scan confirms the phenomenon (a stable, reproducible azimuthal anisotropy), the underlying cause — phase desynchronization — remains an interpretative element of the theoretical framework. Future datasets (CHIME/FRB Catalog 3, DSA-2000, SKA) will be required to test the mechanism more directly.

Why the Columbus Scan Detects an Azimuthal Signal

The anisotropy revealed by the Columbus Rotational Scan is not expected to appear as a classical dipole. In the Siamese Cosmology framework, the predicted imprint is not a surplus of matter or energy along a single direction, but a phase-dependent modulation around a fixed axis.

A phase–desynchronization Δϕ(a) between the two CPT-related cosmic sectors generates a rotational contrast, not a polar gradient. As a consequence, the measurable quantity is the azimuthal variation of DM when rotating a hemispheric divider around the axis. This produces a sinusoidal dependence:
Δ⟨DM⟩(ψ)∝sin(ψ−ϕ0​)

This is exactly the form recovered in the Columbus Scan.

The fact that Mode A (standard dipole test) remains null, while Mode B shows a clean sinusoid with a stable phase ϕ0​ across independent datasets, is precisely the geometric signature expected from this mechanism.

In short: the Columbus-axis modulation is rotational, not dipolar. This explains why its phase is stable and why the sinusoidal pattern persists regardless of catalog heterogeneity.
