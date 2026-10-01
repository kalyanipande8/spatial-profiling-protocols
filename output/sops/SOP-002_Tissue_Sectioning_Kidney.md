# SOP-002: Tissue Sectioning Protocol for Kidney Tissue (CODEX / Xenium Workflows)

**Effective date:** 2026-09-30
**Applies to:** Kidney tissue blocks (FFPE or fresh-frozen/OCT) being
sectioned for CODEX or Xenium multiplexed spatial profiling.
**Derived from:** Protocol optimization study, variant `thickness_3um`
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
| Section thickness | 3.0 um | Highest mean image-quality score across 6 screening replicates (mean = 91.35, sd = 1.92) |
| Measured downstream consistency at this thickness | 91.67% of validation runs >= quality threshold (n=60) | See `output/consistency_summary.csv` |

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
   (3.0 um). Do not deviate more than +/-1 um from this
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
| 1.0 | 2026-09-30 | Initial SOP derived from protocol optimization study |
