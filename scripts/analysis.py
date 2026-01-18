#!/usr/bin/env python3
"""
analysis.py — CHIME/FRB Catalog 2: Columbus Scan + ajuste sinusoidal
para confirmación del eje siamés (RA=170°, Dec=40°).

Estructura esperada de carpetas (relative a este archivo):

chime2/
├── data/
│   └── chimefrbcat2.csv
├── src/
│   └── analysis.py
└── figs/
    └── (salidas PNG, CSV, JSON)

Autor: CosmicThinker & ChatGPT (Toko)
"""

import json
import math
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Barra de progreso: opcional
try:
    from tqdm.auto import tqdm
except ImportError:  # fallback mínimo
    tqdm = None


# =========================
#  CONFIGURACIÓN GLOBAL
# =========================

# Eje siamés (en grados, J2000)
RA_SIAMESE_DEG = 170.0
DEC_SIAMESE_DEG = 40.0

# Filtros estándar
MIN_ABS_B_DEG = 20.0      # |b| > 20°
MIN_DM = 800.0            # DM > 800

# Columbus Scan
SCAN_STEP_DEG = 1.0       # paso en ψ (°), típico 1°
DM_COLUMN_CANDIDATES = ["dm", "dm_excess", "dm_inferred", "dm_exc", "dm_fitb"]

# Permutaciones
N_PERMUTATIONS = 1000     # puedes bajar para tests rápidos
RANDOM_SEED = 12345       # reproducible


# =========================
#  UTILIDADES
# =========================

def get_paths():
    """Devuelve rutas raíz, data y figs basadas en la ubicación de este archivo."""
    this_file = Path(__file__).resolve()
    root = this_file.parent.parent      # .../chime2
    data_dir = root / "data"
    figs_dir = root / "figs"
    figs_dir.mkdir(exist_ok=True)
    return root, data_dir, figs_dir


def find_column(df: pd.DataFrame, candidates, required=True, context=""):
    """
    Busca la primera columna cuyo nombre (case-insensitive) coincide con
    alguno de los candidatos. Si required=True y no encuentra, lanza ValueError.
    """
    cols_lower = {c.lower(): c for c in df.columns}
    for cand in candidates:
        cl = cand.lower()
        if cl in cols_lower:
            return cols_lower[cl]
    if required:
        raise ValueError(
            f"No se encontró ninguna de las columnas {candidates} en el DataFrame "
            f"para {context}. Columnas disponibles: {list(df.columns)}"
        )
    return None


def to_unit_vector(ra_deg, dec_deg):
    """Convierte RA,Dec (deg) a vector unitario cartesiano (x,y,z)."""
    ra_rad = np.deg2rad(ra_deg)
    dec_rad = np.deg2rad(dec_deg)
    cos_dec = np.cos(dec_rad)
    x = cos_dec * np.cos(ra_rad)
    y = cos_dec * np.sin(ra_rad)
    z = np.sin(dec_rad)
    return np.stack((x, y, z), axis=-1)


def build_axis_basis(ra0_deg, dec0_deg):
    """
    Construye una base ortonormal (e1, e2, e3) donde:
      - e3 es la dirección del eje siamés
      - e1, e2 generan el plano perpendicular al eje
    """
    axis = to_unit_vector(ra0_deg, dec0_deg)
    # Asegurar normalización
    axis = axis / np.linalg.norm(axis)

    # Vector temporal no paralelo al eje
    temp = np.array([0.0, 0.0, 1.0])
    if np.allclose(axis, temp, atol=1e-3) or np.allclose(axis, -temp, atol=1e-3):
        temp = np.array([0.0, 1.0, 0.0])

    e1 = np.cross(temp, axis)
    e1 = e1 / np.linalg.norm(e1)
    e2 = np.cross(axis, e1)
    e2 = e2 / np.linalg.norm(e2)
    e3 = axis
    return e1, e2, e3


def compute_phi_deg(ra_deg, dec_deg, ra0_deg=RA_SIAMESE_DEG, dec0_deg=DEC_SIAMESE_DEG):
    """
    Calcula el ángulo φ (0–360°) alrededor del eje siamés para cada (RA,Dec).
    """
    vecs = to_unit_vector(ra_deg, dec_deg)
    e1, e2, e3 = build_axis_basis(ra0_deg, dec0_deg)

    # Proyecciones en la base (e1, e2, e3)
    x_prime = vecs @ e1
    y_prime = vecs @ e2
    # z_prime = vecs @ e3  # no lo usamos para φ

    phi_rad = np.arctan2(y_prime, x_prime)  # rango (-π, π]
    phi_deg = np.rad2deg(phi_rad)
    phi_deg = np.mod(phi_deg, 360.0)  # 0–360
    return phi_deg


# =========================
#  CARGA Y PREPROCESO
# =========================

def load_and_prepare_catalog(csv_path: Path) -> pd.DataFrame:
    """
    Carga el catálogo CHIME/FRB Catalog 2 desde csv_path y devuelve
    un DataFrame con columnas estandarizadas:
      ra_deg, dec_deg, dm, (opcional b_deg) y phi_deg.
    Aplica filtros en |b| y DM si hay datos disponibles.
    """
    print(f"[INFO] Cargando catálogo desde: {csv_path}")
    df = pd.read_csv(csv_path)

    # Detectar columnas de RA y Dec
    ra_col = find_column(df, ["ra", "raj", "ra_deg"], context="RA")
    dec_col = find_column(df, ["dec", "decj", "dec_deg"], context="Dec")

    # Detectar columna DM
    dm_col = find_column(df, DM_COLUMN_CANDIDATES, context="DM")

    # Detectar latitud galáctica (si existe)
    b_col = find_column(df, ["b", "gb", "gal_b_deg"], required=False, context="b (galactic lat)")

    # Renombrar a estándar
    df = df.rename(
        columns={
            ra_col: "ra_deg",
            dec_col: "dec_deg",
            dm_col: "dm",
            **({b_col: "b_deg"} if b_col is not None else {})
        }
    )

    # Eliminar filas con NaN en columnas claves
    required_cols = ["ra_deg", "dec_deg", "dm"]
    before = len(df)
    df = df.dropna(subset=required_cols)
    after = len(df)
    print(f"[INFO] Filas sin NaN en {required_cols}: {after} (eliminadas {before - after})")

    # Filtro en |b| > MIN_ABS_B_DEG si hay b_deg
    if "b_deg" in df.columns and MIN_ABS_B_DEG is not None:
        before = len(df)
        df = df.loc[df["b_deg"].abs() > MIN_ABS_B_DEG].copy()
        after = len(df)
        print(f"[INFO] Filtro |b| > {MIN_ABS_B_DEG}°: quedan {after} eventos (eliminadas {before - after})")
    else:
        print("[WARN] No se encontró b_deg; se omite el filtro en latitud galáctica.")

    # Filtro DM > MIN_DM
    if MIN_DM is not None:
        before = len(df)
        df = df.loc[df["dm"] > MIN_DM].copy()
        after = len(df)
        print(f"[INFO] Filtro DM > {MIN_DM}: quedan {after} eventos (eliminadas {before - after})")

    # Calcular φ respecto al eje siamés
    print("[INFO] Calculando φ (phi_deg) respecto al eje siamés...")
    df["phi_deg"] = compute_phi_deg(df["ra_deg"].values, df["dec_deg"].values)

    print(f"[INFO] Catálogo preparado con {len(df)} FRBs.")
    return df


# =========================
#  COLUMBUS SCAN (MODO B)
# =========================

def columbus_scan_mode_B(df: pd.DataFrame, step_deg: float = SCAN_STEP_DEG) -> pd.DataFrame:
    """
    Columbus Scan modo B:
      - Se recorre ψ en [0, 360) con paso step_deg.
      - Para cada ψ:
           H1: φ ∈ [ψ, ψ+180)
           H2: φ ∈ [ψ+180, ψ+360)
        (mod 360)
      - Se calcula Δ⟨DM⟩(ψ) = mean(dm_H1) - mean(dm_H2).
    Devuelve DataFrame con columnas:
      psi_deg, delta_dm, n_hemi1, n_hemi2
    """
    if "phi_deg" not in df.columns:
        raise ValueError("La columna 'phi_deg' no está en el DataFrame. Ejecuta primero el cálculo de φ.")

    phi = df["phi_deg"].values
    dm = df["dm"].values

    psi_values = np.arange(0.0, 360.0, step_deg)
    delta_list = []
    n1_list = []
    n2_list = []

    print(f"[INFO] Ejecutando Columbus Scan modo B con paso {step_deg}°...")

    for psi in psi_values:
        # H1: φ en [psi, psi+180) modulo 360
        # Convertimos a rango 0–360
        phi_shift = (phi - psi) % 360.0
        mask_h1 = (phi_shift >= 0.0) & (phi_shift < 180.0)
        mask_h2 = ~mask_h1

        dm_h1 = dm[mask_h1]
        dm_h2 = dm[mask_h2]

        if len(dm_h1) == 0 or len(dm_h2) == 0:
            delta = np.nan
        else:
            delta = dm_h1.mean() - dm_h2.mean()

        delta_list.append(delta)
        n1_list.append(len(dm_h1))
        n2_list.append(len(dm_h2))

    scan_df = pd.DataFrame(
        {
            "psi_deg": psi_values,
            "delta_dm": delta_list,
            "n_hemi1": n1_list,
            "n_hemi2": n2_list,
        }
    )

    # Eliminar posibles NaN si hay ángulos sin datos en un hemisferio
    before = len(scan_df)
    scan_df = scan_df.dropna(subset=["delta_dm"])
    after = len(scan_df)
    if after < before:
        print(f"[WARN] Eliminados {before - after} puntos del scan con delta_dm = NaN (hemisferios vacíos).")

    print(f"[INFO] Columbus Scan completado con {len(scan_df)} puntos válidos.")
    return scan_df


# =========================
#  AJUSTE SINUSOIDAL (LINEAL)
# =========================

def fit_sine_linear(scan_df: pd.DataFrame):
    """
    Ajuste lineal de:
        y = B sin(phi) + D cos(phi) + C
    donde phi = psi_deg (en radianes) y y = delta_dm.

    Devuelve:
      - params: dict con B, D, C, A_eff, phi0_deg, R2
    """
    phi_deg = scan_df["psi_deg"].values
    y = scan_df["delta_dm"].values

    phi_rad = np.deg2rad(phi_deg)
    X = np.column_stack([np.sin(phi_rad), np.cos(phi_rad), np.ones_like(phi_rad)])

    # Ajuste lineal por mínimos cuadrados
    # Resolviendo X beta = y
    beta, residuals, rank, s = np.linalg.lstsq(X, y, rcond=None)
    B, D, C = beta

    # Predicción y R^2
    y_pred = X @ beta
    ss_res = np.sum((y - y_pred) ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    R2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan

    # Convertir (B, D) en amplitud y fase:
    # A_eff sin(phi - phi0) = B sin(phi) + D cos(phi)
    # => A_eff = sqrt(B^2 + D^2)
    #    phi0 = atan2(D, B)
    A_eff = math.sqrt(B**2 + D**2)
    phi0_rad = math.atan2(D, B)
    phi0_deg = math.degrees(phi0_rad)  # puede estar en (-180,180]
    # Normalizar a [0, 360)
    phi0_deg = phi0_deg % 360.0

    params = {
        "B": float(B),
        "D": float(D),
        "C": float(C),
        "A_eff": float(A_eff),
        "phi0_deg": float(phi0_deg),
        "R2": float(R2),
    }
    return params, y_pred


# =========================
#  PERMUTATION TEST
# =========================

def permutation_test_amplitude(df: pd.DataFrame,
                               step_deg: float,
                               n_permutations: int = N_PERMUTATIONS,
                               random_seed: int = RANDOM_SEED):
    """
    Test de permutación sobre la amplitud |A_eff| del ajuste senoidal
    del Columbus Scan.

    1) Calcula |A_eff|_data a partir del df original.
    2) Para cada permutación:
         - baraja los valores de DM entre FRBs
         - recalcula scan + ajuste
         - almacena |A_eff|_perm
    3) Devuelve:
         - p_perm = frac(|A_perm| >= |A_data|)
         - lista de amplitudes permutadas

    Usa tqdm si está disponible para barra de progreso.
    """
    rng = np.random.default_rng(random_seed)

    # Amplitud real
    scan_df = columbus_scan_mode_B(df, step_deg=step_deg)
    params_data, _ = fit_sine_linear(scan_df)
    A_data = abs(params_data["A_eff"])
    print(f"[INFO] Amplitud real |A_eff| = {A_data:.3f} pc cm^-3")

    amps_perm = []

    iterator = range(n_permutations)
    if tqdm is not None:
        iterator = tqdm(iterator, desc="Permutaciones", unit="perm")

    dm_original = df["dm"].values.copy()

    for _ in iterator:
        dm_shuffled = dm_original.copy()
        rng.shuffle(dm_shuffled)
        df_perm = df.copy()
        df_perm["dm"] = dm_shuffled

        scan_perm = columbus_scan_mode_B(df_perm, step_deg=step_deg)
        params_perm, _ = fit_sine_linear(scan_perm)
        A_perm = abs(params_perm["A_eff"])
        amps_perm.append(A_perm)

    amps_perm = np.array(amps_perm)
    p_perm = float(np.mean(amps_perm >= A_data))

    print(f"[INFO] Test de permutación completado: p_perm(|A_eff|) = {p_perm:.5f}")
    return A_data, amps_perm, p_perm


# =========================
#  PLOTS
# =========================

def plot_phi_histogram(df: pd.DataFrame, outpath: Path):
    plt.figure(figsize=(8, 4))
    plt.hist(df["phi_deg"].values, bins=36)
    plt.xlabel("φ (deg) respecto al eje siamés")
    plt.ylabel("Número de FRBs")
    plt.title("Distribución de φ — CHIME/FRB Catalog 2")
    plt.tight_layout()
    plt.savefig(outpath, dpi=150)
    plt.close()
    print(f"[INFO] Histograma de φ guardado en: {outpath}")


def plot_columbus_scan(scan_df: pd.DataFrame, outpath: Path):
    plt.figure(figsize=(8, 4))
    plt.plot(scan_df["psi_deg"].values, scan_df["delta_dm"].values, marker=".", linestyle="-")
    plt.xlabel("ψ (deg) — ángulo de corte hemisférico")
    plt.ylabel("Δ⟨DM⟩ (pc cm$^{-3}$)")
    plt.title("Columbus Scan (modo B) — CHIME/FRB Catalog 2")
    plt.tight_layout()
    plt.savefig(outpath, dpi=150)
    plt.close()
    print(f"[INFO] Columbus Scan guardado en: {outpath}")


def plot_columbus_with_fit(scan_df: pd.DataFrame, y_fit, params: dict, outpath: Path):
    plt.figure(figsize=(8, 4))
    psi = scan_df["psi_deg"].values
    plt.plot(psi, scan_df["delta_dm"].values, marker=".", linestyle="", label="Datos")
    plt.plot(psi, y_fit, linestyle="-", label="Ajuste senoidal")
    plt.xlabel("ψ (deg)")
    plt.ylabel("Δ⟨DM⟩ (pc cm$^{-3}$)")
    title = (
        f"Columbus Scan + sin-fit — A={params['A_eff']:.1f}, "
        f"φ₀={params['phi0_deg']:.1f}°, R²={params['R2']:.2f}"
    )
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.savefig(outpath, dpi=150)
    plt.close()
    print(f"[INFO] Columbus Scan + ajuste guardado en: {outpath}")


# =========================
#  MAIN
# =========================

def main():
    root, data_dir, figs_dir = get_paths()
    csv_path = data_dir / "chimefrbcat2.csv"

    if not csv_path.exists():
        print(f"[ERROR] No se encontró el archivo: {csv_path}")
        sys.exit(1)

    # 1) Cargar y preparar catálogo
    df = load_and_prepare_catalog(csv_path)

    # 2) Histograma de φ
    plot_phi_histogram(df, figs_dir / "phi_distribution.png")

    # 3) Columbus Scan
    scan_df = columbus_scan_mode_B(df, step_deg=SCAN_STEP_DEG)

    # Guardar CSV con el scan
    scan_csv = figs_dir / "chime2_columbus_scan.csv"
    scan_df.to_csv(scan_csv, index=False)
    print(f"[INFO] Columbus Scan CSV guardado en: {scan_csv}")

    # 4) Ajuste senoidal (lineal)
    params, y_fit = fit_sine_linear(scan_df)
    print("[INFO] Parámetros del ajuste senoidal (lineal):")
    for k, v in params.items():
        print(f"   {k}: {v}")

    # 5) Gráficas del scan
    plot_columbus_scan(scan_df, figs_dir / "columbus_scan.png")
    plot_columbus_with_fit(scan_df, y_fit, params, figs_dir / "columbus_scan_with_fit.png")

    # 6) Test de permutación
    A_data, amps_perm, p_perm = permutation_test_amplitude(
        df,
        step_deg=SCAN_STEP_DEG,
        n_permutations=N_PERMUTATIONS,
        random_seed=RANDOM_SEED,
    )

    # 7) Guardar resumen en JSON
    summary = {
        "timestamp_utc": datetime.utcnow().isoformat() + "Z",
        "catalog": "CHIME/FRB Catalog 2",
        "input_csv": str(csv_path),
        "filters": {
            "min_abs_b_deg": MIN_ABS_B_DEG,
            "min_dm": MIN_DM,
        },
        "siamese_axis": {
            "ra_deg": RA_SIAMESE_DEG,
            "dec_deg": DEC_SIAMESE_DEG,
        },
        "scan": {
            "step_deg": SCAN_STEP_DEG,
            "n_points": int(len(scan_df)),
        },
        "fit_params": params,
        "permutation_test": {
            "n_permutations": N_PERMUTATIONS,
            "A_data_abs": float(A_data),
            "p_perm_absA": float(p_perm),
        },
    }

    summary_path = figs_dir / "chime2_sinefit_summary.json"
    with open(summary_path, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)
    print(f"[INFO] Resumen del ajuste + permutaciones guardado en: {summary_path}")

    # Guardar amplitudes permutadas (por si queremos histograma después)
    amps_perm_path = figs_dir / "chime2_perm_amplitudes.npy"
    np.save(amps_perm_path, amps_perm)
    print(f"[INFO] Amplitudes permutadas guardadas en: {amps_perm_path}")

    print("\n[DONE] Pipeline completo ejecutado.")
    print("Resultados clave en 'figs/' (PNG, CSV, JSON).")


if __name__ == "__main__":
    main()
