"""
generate_sops.py

Generates two Standard Operating Procedure documents under output/sops/,
populated with the real optimized protocol parameters and the real
statistics produced by analysis.py:
  - SOP-001_HE_Staining_Kidney.md
  - SOP-002_Tissue_Sectioning_Kidney.md

Must be run after analysis.py (it reads output/variant_comparison.csv and
output/consistency_summary.csv).
"""

import os
from datetime import date

import pandas as pd

PROJECT_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.path.join(PROJECT_DIR, "output")
SOP_DIR = os.path.join(OUTPUT_DIR, "sops")


def _load_optimal_params():
    ranked = pd.read_csv(os.path.join(OUTPUT_DIR, "variant_comparison.csv"))
    top = ranked.iloc[0]
    consistency = pd.read_csv(os.path.join(OUTPUT_DIR, "consistency_summary.csv"))
    overall = consistency.loc[consistency["platform"] == "overall"].iloc[0]
    return top, overall, consistency


def write_he_staining_sop(top, overall):
    text = f"""# SOP-001: H&E Staining Protocol for Kidney Tissue Sections (CODEX / Xenium Workflows)

**Effective date:** {date.today().isoformat()}
**Applies to:** Formalin-fixed paraffin-embedded (FFPE) and fresh-frozen kidney
tissue sections prior to CODEX or Xenium multiplexed imaging.
**Derived from:** Protocol optimization study, variant `{top['variant_id']}`
(see `output/analysis_report.md` and `output/variant_comparison.csv`).

> Synthetic-data disclaimer: the optimized parameters below were derived from
> a simulated experiment set (`simulate_protocol_experiments.py`), not from
> real patient tissue or a real production assay. The document format and
> the optimization methodology are real and reusable; the specific numbers
> are illustrative.

## 1. Purpose

Define the staining-time and antibody-concentration parameters that
maximize downstream image-quality consistency (signal-to-noise, staining
uniformity, tissue integrity, marker-channel alignment) for kidney tissue
processed through CODEX or Xenium multiplexed imaging.

## 2. Optimized Parameters

| Parameter | Optimized value | Basis |
|---|---|---|
| Primary/probe staining incubation time | {int(top['staining_time_min'])} min | Highest mean image-quality score across {int(top['n'])} screening replicates (mean = {top['mean_quality']:.2f}, sd = {top['std_quality']:.2f}) |
| Antibody / probe concentration | {top['antibody_conc_ug_ml']:.1f} ug/mL | Same screening comparison |
| Downstream validation run size | {int(overall['n'])} replicate experiments across CODEX and Xenium | Consistency confirmation run |
| Measured image-quality consistency at these settings | {overall['pct_consistent']:.2f}% of runs >= quality threshold | See `output/consistency_summary.csv` |

## 3. Materials

- Xylene / clearing agent (for FFPE) or fresh-frozen OCT-embedded kidney block
- Ethanol series (100%, 95%, 70%)
- Hematoxylin (Mayer's or Harris)
- Eosin Y, alcoholic
- Acid alcohol differentiation solution (0.3% HCl in 70% ethanol)
- Bluing reagent (ammonia water or equivalent)
- Mounting medium and coverslips
- Humidity-controlled staining chamber

## 4. Procedure

1. Section kidney tissue per SOP-002 (Tissue Sectioning) at the validated
   thickness before proceeding.
2. Deparaffinize/rehydrate (FFPE) or fix fresh-frozen sections in cold
   acetone/ethanol per standard rehydration steps.
3. Stain in hematoxylin for the optimized incubation window
   ({int(top['staining_time_min'])} min at room temperature); staining time outside
   this window measurably reduces staining uniformity and SNR in the
   optimization data (see Section 5).
4. Rinse in running deionized water until clear.
5. Differentiate in acid alcohol (1-2 quick dips), then blue.
6. Counterstain in eosin Y per standard timing for the tissue type.
7. Dehydrate through the ethanol series, clear, and mount.
8. Proceed immediately to antibody/probe hybridization at the optimized
   concentration ({top['antibody_conc_ug_ml']:.1f} ug/mL) for CODEX antibody panels or
   Xenium probe panels.
9. Image and score the run against the four image-quality metrics defined
   in `simulate_protocol_experiments.py` (SNR, staining uniformity, tissue
   integrity, marker-channel alignment) before releasing data for analysis.

## 5. Rationale (from optimization data)

- Staining incubation time has a **quadratic** effect on both SNR and
  staining uniformity in the screening data: deviating from the optimum in
  either direction (under- or over-staining) reduces both metrics.
  Over-staining beyond ~10 minutes past optimum additionally reduces tissue
  integrity.
- Antibody/probe concentration also has a quadratic effect: too little
  under-labels (low SNR), too much raises non-specific background and
  patchiness (lower uniformity).
- A Kruskal-Wallis test across all screening variants and a Welch's
  t-test of the best vs. worst variant (both reported in
  `output/analysis_report.md`) confirm the optimized variant significantly
  outperforms poorly-parameterized alternatives.

## 6. Acceptance Criteria

A run is accepted for downstream analysis if its composite image-quality
score is >= the validated spec threshold used in `analysis.py`
(`QUALITY_PASS_THRESHOLD`). Runs failing this threshold should be
re-sectioned and re-stained rather than analyzed.

## 7. Revision History

| Version | Date | Change |
|---|---|---|
| 1.0 | {date.today().isoformat()} | Initial SOP derived from protocol optimization study |
"""
    os.makedirs(SOP_DIR, exist_ok=True)
    path = os.path.join(SOP_DIR, "SOP-001_HE_Staining_Kidney.md")
    with open(path, "w") as f:
        f.write(text)
    return path


def write_sectioning_sop(top, overall):
    text = f"""# SOP-002: Tissue Sectioning Protocol for Kidney Tissue (CODEX / Xenium Workflows)

**Effective date:** {date.today().isoformat()}
**Applies to:** Kidney tissue blocks (FFPE or fresh-frozen/OCT) being
sectioned for CODEX or Xenium multiplexed spatial profiling.
**Derived from:** Protocol optimization study, variant `{top['variant_id']}`
(see `output/analysis_report.md` and `output/variant_comparison.csv`).

> Synthetic-data disclaimer: the optimized parameters below were derived from
> a simulated experiment set (`simulate_protocol_experiments.py`), not from
> real patient tissue or a real production assay. The document format and
> the optimization methodology are real and reusable; the specific numbers
> are illustrative.

## 1. Purpose

Define the sectioning thickness that maximizes tissue-integrity and
marker/probe-channel alignment for kidney tissue carried through
multi-cycle CODEX imaging or single-round Xenium in-situ imaging.

## 2. Optimized Parameter

| Parameter | Optimized value | Basis |
|---|---|---|
| Section thickness | {top['thickness_um']:.1f} um | Highest mean image-quality score across {int(top['n'])} screening replicates (mean = {top['mean_quality']:.2f}, sd = {top['std_quality']:.2f}) |
| Measured downstream consistency at this thickness | {overall['pct_consistent']:.2f}% of validation runs >= quality threshold (n={int(overall['n'])}) | See `output/consistency_summary.csv` |

## 3. Materials

- Calibrated rotary microtome (FFPE) or cryostat (fresh-frozen/OCT)
- Positively charged glass slides
- Water bath / ribbon-flotation bath (FFPE) set to standard temperature
- Anti-roll plate (cryostat sectioning)
- Slide warmer / drying rack

## 4. Procedure

1. Equilibrate the tissue block to the microtome/cryostat chamber
   temperature per standard practice for kidney tissue.
2. Trim the block face until a full, even tissue profile is exposed.
3. Cut serial sections at the optimized thickness
   ({top['thickness_um']:.1f} um). Do not deviate more than +/-1 um from this
   target - the optimization data show tissue-integrity and
   marker-channel-alignment scores degrade roughly linearly with distance
   from this thickness, and the effect is more pronounced for CODEX
   (multi-cycle) than Xenium (single-round) imaging.
4. Float sections onto the water bath (FFPE) or pick up directly with a
   charged slide (cryostat), avoiding folds, tears, or compression.
5. Mount on positively charged slides and dry flat.
6. Proceed to H&E staining (SOP-001) or directly to antibody/probe
   application, depending on workflow.
7. Visually QC each section for folds, tears, or incomplete tissue
   transfer before proceeding; reject and re-cut sections that fail QC.

## 5. Rationale (from optimization data)

- Sectioning thickness deviation (either thinner or thicker than optimum)
  reduces tissue-integrity score roughly linearly with |deviation|.
- Sections thicker than optimum specifically reduce marker/probe-channel
  alignment, because greater z-depth increases drift across imaging
  cycles; this effect is larger for CODEX (repeated antibody
  stripping/re-staining cycles) than for Xenium (single imaging round).
- The optimized thickness was selected via a Kruskal-Wallis comparison
  across all screening variants (thickness, staining time, and antibody
  concentration swept independently) and confirmed by a Welch's t-test of
  the best vs. worst-performing variant - both reported in
  `output/analysis_report.md`.

## 6. Acceptance Criteria

A section is accepted for staining/imaging if it meets visual QC (no
folds/tears, full tissue transfer) and the resulting image-quality score
is >= the validated spec threshold used in `analysis.py`
(`QUALITY_PASS_THRESHOLD`).

## 7. Revision History

| Version | Date | Change |
|---|---|---|
| 1.0 | {date.today().isoformat()} | Initial SOP derived from protocol optimization study |
"""
    os.makedirs(SOP_DIR, exist_ok=True)
    path = os.path.join(SOP_DIR, "SOP-002_Tissue_Sectioning_Kidney.md")
    with open(path, "w") as f:
        f.write(text)
    return path


def main():
    top, overall, consistency = _load_optimal_params()
    he_path = write_he_staining_sop(top, overall)
    sectioning_path = write_sectioning_sop(top, overall)
    print(f"Wrote {he_path}")
    print(f"Wrote {sectioning_path}")
    return he_path, sectioning_path


if __name__ == "__main__":
    main()
