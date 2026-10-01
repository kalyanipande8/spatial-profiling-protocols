"""
simulate_protocol_experiments.py

Simulates two stages of CODEX / Xenium multiplexed tissue-imaging
experiments on kidney tissue, each producing real computed image-quality
metrics from an underlying (synthetic) staining/sectioning physics model
plus measurement noise. Nothing here is hand-picked after the fact - the
"optimal" protocol is only known to this script as the center of the
simulated response surface; `analysis.py` has to rediscover it from the
data, the same way a real protocol-optimization study would.

Stage 1 - SCREENING (one-factor-at-a-time design):
  A baseline protocol (30 min staining, 5 ug/mL antibody concentration,
  5 um sectioning thickness) is run on both platforms, plus 3 additional
  levels for each of the 3 factors (holding the other two at baseline).
  10 protocol variants x 2 platforms x 3 replicates = 60 experiments.

Stage 2 - VALIDATION:
  Once `analysis.py` identifies the best-performing variant from Stage 1,
  that single protocol is re-run 30 times on each platform (60 experiments
  total) to characterize its run-to-run image-quality consistency.

Four image-quality metrics are computed per experiment, each with its own
noise-corrupted dependence on the protocol parameters and platform:
  - snr_db                    signal-to-noise ratio (dB)
  - staining_uniformity_pct   coefficient-of-variation-based uniformity (0-100)
  - tissue_integrity_pct      morphology/RNA-or-epitope-integrity score (0-100)
  - marker_alignment_pct      multi-cycle / multi-probe channel co-registration (0-100)

A single composite `image_quality_score` (0-100) is the weighted average of
the four metrics (SNR min-max scaled to 0-100 first). A fixed random seed
makes every run reproducible.
"""

import os

import numpy as np
import pandas as pd

PROJECT_DIR = os.path.dirname(os.path.abspath(__file__))
RAW_DIR = os.path.join(PROJECT_DIR, "data", "raw")

RANDOM_SEED = 42
PLATFORMS = ["CODEX", "Xenium"]

# "True" (unknown-to-the-analysis) optimum of the simulated response surface.
TRUE_OPTIMUM = dict(staining_time_min=30, antibody_conc_ug_ml=5.0, thickness_um=5.0)

BASELINE = dict(staining_time_min=30, antibody_conc_ug_ml=5.0, thickness_um=5.0)

# One-factor-at-a-time screening levels (baseline + 3 alternate levels per factor).
SCREENING_VARIANTS = [
    dict(variant_id="baseline", staining_time_min=30, antibody_conc_ug_ml=5.0, thickness_um=5.0),
    # staining time sweep (antibody conc / thickness held at baseline)
    dict(variant_id="time_15min", staining_time_min=15, antibody_conc_ug_ml=5.0, thickness_um=5.0),
    dict(variant_id="time_45min", staining_time_min=45, antibody_conc_ug_ml=5.0, thickness_um=5.0),
    dict(variant_id="time_60min", staining_time_min=60, antibody_conc_ug_ml=5.0, thickness_um=5.0),
    # antibody concentration sweep
    dict(variant_id="conc_2ugml", staining_time_min=30, antibody_conc_ug_ml=2.0, thickness_um=5.0),
    dict(variant_id="conc_10ugml", staining_time_min=30, antibody_conc_ug_ml=10.0, thickness_um=5.0),
    dict(variant_id="conc_20ugml", staining_time_min=30, antibody_conc_ug_ml=20.0, thickness_um=5.0),
    # sectioning thickness sweep
    dict(variant_id="thickness_3um", staining_time_min=30, antibody_conc_ug_ml=5.0, thickness_um=3.0),
    dict(variant_id="thickness_7um", staining_time_min=30, antibody_conc_ug_ml=5.0, thickness_um=7.0),
    dict(variant_id="thickness_10um", staining_time_min=30, antibody_conc_ug_ml=5.0, thickness_um=10.0),
]

N_REPLICATES_SCREENING = 3
N_REPLICATES_VALIDATION = 30  # per platform -> 60 total validation experiments

SNR_MIN, SNR_MAX = 8.0, 32.0  # min-max scaling range used to fold SNR into the composite score
QUALITY_WEIGHTS = dict(snr=0.30, uniformity=0.25, integrity=0.25, alignment=0.20)
QUALITY_PASS_THRESHOLD = 87.5  # composite image_quality_score >= this counts as "in spec"


def _compute_metrics(rng, staining_time_min, antibody_conc_ug_ml, thickness_um, platform):
    """Real (simulated) image-quality metrics from a synthetic physics model.

    Deviations from TRUE_OPTIMUM degrade each metric in a mechanistically
    plausible way, plus platform-specific offsets and per-experiment noise.
    """
    time_dev = staining_time_min - TRUE_OPTIMUM["staining_time_min"]
    conc_dev = antibody_conc_ug_ml - TRUE_OPTIMUM["antibody_conc_ug_ml"]
    thick_dev = thickness_um - TRUE_OPTIMUM["thickness_um"]
    platform_codex = 1.0 if platform == "CODEX" else 0.0

    # SNR (dB): too little antibody starves the signal, too much raises background;
    # under-staining also hurts; Xenium's probe-based chemistry runs slightly hotter SNR.
    snr_mean = (
        27.5
        - 0.055 * conc_dev ** 2
        - 0.020 * time_dev ** 2
        - 0.35 * platform_codex
    )
    snr = float(np.clip(rng.normal(snr_mean, 1.3), 1.0, None))

    # Staining uniformity (%): sensitive to both under- and over-staining time and to
    # antibody concentration drifting off-optimum (patchy at low conc, saturating/uneven at high).
    uniformity_mean = (
        96.5
        - 0.065 * time_dev ** 2
        - 0.09 * conc_dev ** 2
    )
    uniformity = float(np.clip(rng.normal(uniformity_mean, 1.8), 0.0, 100.0))

    # Tissue integrity (%): thinner/thicker-than-optimal sections tear or fold, and
    # prolonged staining beyond optimum degrades morphology/epitopes/RNA.
    integrity_mean = (
        97.0
        - 0.95 * abs(thick_dev)
        - 0.10 * max(0.0, time_dev - 10)
    )
    integrity = float(np.clip(rng.normal(integrity_mean, 1.6), 0.0, 100.0))

    # Marker-channel alignment (%): thicker sections cause more z-drift across imaging
    # cycles/rounds; CODEX (multi-cycle antibody stripping) is more sensitive to this
    # than Xenium (single-round in-situ hybridization imaging).
    alignment_mean = (
        95.5
        - (1.15 + 0.5 * platform_codex) * max(0.0, thick_dev)
        - 0.6 * max(0.0, -thick_dev)
        - 1.0 * platform_codex
    )
    alignment = float(np.clip(rng.normal(alignment_mean, 1.7), 0.0, 100.0))

    snr_scaled = float(np.clip((snr - SNR_MIN) / (SNR_MAX - SNR_MIN) * 100.0, 0.0, 100.0))
    composite = (
        QUALITY_WEIGHTS["snr"] * snr_scaled
        + QUALITY_WEIGHTS["uniformity"] * uniformity
        + QUALITY_WEIGHTS["integrity"] * integrity
        + QUALITY_WEIGHTS["alignment"] * alignment
    )

    return dict(
        snr_db=round(snr, 2),
        staining_uniformity_pct=round(uniformity, 2),
        tissue_integrity_pct=round(integrity, 2),
        marker_alignment_pct=round(alignment, 2),
        image_quality_score=round(composite, 2),
    )


def simulate_screening():
    rng = np.random.default_rng(RANDOM_SEED)
    rows = []
    exp_id = 1
    for variant in SCREENING_VARIANTS:
        for platform in PLATFORMS:
            for rep in range(1, N_REPLICATES_SCREENING + 1):
                metrics = _compute_metrics(
                    rng, variant["staining_time_min"], variant["antibody_conc_ug_ml"],
                    variant["thickness_um"], platform,
                )
                rows.append({
                    "experiment_id": f"SCR-{exp_id:04d}",
                    "stage": "screening",
                    "variant_id": variant["variant_id"],
                    "platform": platform,
                    "replicate": rep,
                    "tissue_type": "kidney",
                    "staining_time_min": variant["staining_time_min"],
                    "antibody_conc_ug_ml": variant["antibody_conc_ug_ml"],
                    "thickness_um": variant["thickness_um"],
                    **metrics,
                })
                exp_id += 1

    df = pd.DataFrame(rows)
    os.makedirs(RAW_DIR, exist_ok=True)
    df.to_csv(os.path.join(RAW_DIR, "screening_experiments.csv"), index=False)
    print(f"Simulated {len(df)} screening experiments "
          f"({len(SCREENING_VARIANTS)} variants x {len(PLATFORMS)} platforms x {N_REPLICATES_SCREENING} reps)")
    return df


def simulate_validation(optimal_params):
    """Re-run the optimized protocol many times per platform to measure consistency.

    `optimal_params` must have keys staining_time_min, antibody_conc_ug_ml,
    thickness_um - normally the output of analysis.find_optimal_variant().
    """
    # Offset the seed from screening so validation draws are independent, but still fixed.
    rng = np.random.default_rng(RANDOM_SEED + 1000)
    rows = []
    exp_id = 1
    for platform in PLATFORMS:
        for rep in range(1, N_REPLICATES_VALIDATION + 1):
            metrics = _compute_metrics(
                rng, optimal_params["staining_time_min"], optimal_params["antibody_conc_ug_ml"],
                optimal_params["thickness_um"], platform,
            )
            rows.append({
                "experiment_id": f"VAL-{exp_id:04d}",
                "stage": "validation",
                "variant_id": "optimized",
                "platform": platform,
                "replicate": rep,
                "tissue_type": "kidney",
                "staining_time_min": optimal_params["staining_time_min"],
                "antibody_conc_ug_ml": optimal_params["antibody_conc_ug_ml"],
                "thickness_um": optimal_params["thickness_um"],
                **metrics,
            })
            exp_id += 1

    df = pd.DataFrame(rows)
    os.makedirs(RAW_DIR, exist_ok=True)
    df.to_csv(os.path.join(RAW_DIR, "validation_experiments.csv"), index=False)
    print(f"Simulated {len(df)} validation experiments on the optimized protocol "
          f"({len(PLATFORMS)} platforms x {N_REPLICATES_VALIDATION} reps)")
    return df


if __name__ == "__main__":
    simulate_screening()
    simulate_validation(BASELINE)
