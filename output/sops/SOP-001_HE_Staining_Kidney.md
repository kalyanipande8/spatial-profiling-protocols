# SOP-001: H&E Staining Protocol for Kidney Tissue Sections (CODEX / Xenium Workflows)

**Effective date:** 2026-09-30
**Applies to:** Formalin-fixed paraffin-embedded (FFPE) and fresh-frozen kidney
tissue sections prior to CODEX or Xenium multiplexed imaging.
**Derived from:** Protocol optimization study, variant `thickness_3um`
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
| Primary/probe staining incubation time | 30 min | Highest mean image-quality score across 6 screening replicates (mean = 91.35, sd = 1.92) |
| Antibody / probe concentration | 5.0 ug/mL | Same screening comparison |
| Downstream validation run size | 60 replicate experiments across CODEX and Xenium | Consistency confirmation run |
| Measured image-quality consistency at these settings | 91.67% of runs >= quality threshold | See `output/consistency_summary.csv` |

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
   (30 min at room temperature); staining time outside
   this window measurably reduces staining uniformity and SNR in the
   optimization data (see Section 5).
4. Rinse in running deionized water until clear.
5. Differentiate in acid alcohol (1-2 quick dips), then blue.
6. Counterstain in eosin Y per standard timing for the tissue type.
7. Dehydrate through the ethanol series, clear, and mount.
8. Proceed immediately to antibody/probe hybridization at the optimized
   concentration (5.0 ug/mL) for CODEX antibody panels or
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
| 1.0 | 2026-09-30 | Initial SOP derived from protocol optimization study |
