# Spatial Profiling Protocols: CODEX/Xenium Kidney Image-Quality Optimization

A **protocol-optimization and SOP-development** pipeline for multiplexed
tissue-imaging (CODEX / Xenium) on kidney tissue: it screens staining time,
antibody/probe concentration, and sectioning-thickness variants, picks the
best-performing protocol by real statistical comparison, measures its
run-to-run image-quality consistency, and generates the resulting Standard
Operating Procedures.

> **Not to be confused with** the sibling project `spatial-profiling-analysis`
> (segmentation/cell-classification pipelines on spatial imaging data). This
> project is upstream of that: it is about **wet-lab protocol optimization**
> (staining/sectioning parameters, SOPs, image-quality consistency), not
> about processing already-acquired images into cell calls.

> **Disclaimer:** This is a portfolio/demo project. All experiment data is
> **synthetically simulated** by `simulate_protocol_experiments.py` from a
> parametric, literature-plausible response-surface model (staining time,
> antibody/probe concentration, and sectioning thickness each degrade SNR,
> staining uniformity, tissue integrity, and marker-channel alignment in a
> mechanistically motivated way, plus noise). No real patient tissue, real
> instrument output, or real employer data is used. The optimization
> methodology, statistics, and SOP structure are real and reusable.

## What it does

1. **`simulate_protocol_experiments.py`** - Stage 1 (screening): simulates
   a one-factor-at-a-time design on kidney tissue - a baseline protocol
   (30 min staining, 5 ug/mL antibody/probe, 5 um sections) plus 3
   alternate levels for each of the 3 factors, each run on both CODEX and
   Xenium with 3 replicates (**60 screening experiments**). Each experiment
   gets four independently-modeled, noise-corrupted image-quality metrics:
   signal-to-noise ratio (SNR, dB), staining uniformity (%), tissue
   integrity (%), and marker/probe-channel alignment (%), folded into one
   composite `image_quality_score` (0-100).

2. **`analysis.find_optimal_variant()`** picks the screening variant with
   the highest mean `image_quality_score`.

3. **`simulate_protocol_experiments.simulate_validation()`** - Stage 2:
   re-runs that single optimized protocol **30 times on each of CODEX and
   Xenium (60 validation experiments)** to characterize real run-to-run
   consistency.

4. **`standardize.py`** concatenates screening + validation into one
   unified schema (`data/processed/unified_experiments.csv`) and flags
   each experiment against a fixed quality spec.

5. **`analysis.py`** runs a Kruskal-Wallis test across all screening
   variants, a Welch's t-test of the best vs. worst variant, and computes
   the **image-quality-consistency metric**: the percentage of validation
   experiments (optimized protocol) meeting the quality-score spec.

6. **`generate_sops.py`** writes two real SOP documents to
   `output/sops/`, populated with the optimized parameters and the actual
   statistics above:
   - `SOP-001_HE_Staining_Kidney.md`
   - `SOP-002_Tissue_Sectioning_Kidney.md`

7. **`run_all.py`** runs the whole pipeline end-to-end.

## How to run

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
python run_all.py
pytest -q
```

## Results

_All numbers below are copied directly from an actual run of `python
run_all.py` - nothing is hardcoded. Fixed random seeds make them
reproducible._

- **120 total experiments** (60 screening + 60 validation) across CODEX
  and Xenium on kidney tissue.
- **Kruskal-Wallis across 10 screening variants: H = 49.026, p = 1.64e-07**
  (n=60) - protocol parameters have a highly significant effect on image
  quality.
- **Optimal protocol identified: `thickness_3um`** - 30 min staining,
  5.0 ug/mL antibody/probe concentration, 3.0 um sectioning thickness
  (mean image_quality_score = 91.35, sd = 1.92, n=6 screening replicates).
  This essentially ties the original baseline (5 um sections, mean =
  91.18) - both sit near the true optimum of the simulated response
  surface, which is the expected outcome of a well-designed screen.
- **Worst protocol: `time_60min`** (60 min staining; mean = 53.28, sd =
  1.18) - prolonged staining time overwhelmingly dominates the other
  factors in degrading image quality in this model.
- **Best vs. worst (Welch's t-test): t = 41.337, p = 6.40e-11** - the
  optimized protocol's advantage is not noise; it is a >38-point mean
  difference in composite image-quality score.
- **Image-quality consistency on the optimized protocol: 91.67%** of the
  60 validation experiments met the quality-score spec (>= 87.5 /100).
  - CODEX: 90.00% (27/30 in spec, mean score 90.37, sd 1.66)
  - Xenium: 93.33% (28/30 in spec, mean score 91.17, sd 1.67)
  - This lands just under the ~95% target reported in the role summary;
    reported honestly as computed rather than rounded up to hit the
    target exactly.

(Exact values in `output/variant_comparison.csv` and
`output/consistency_summary.csv`; full narrative in
`output/analysis_report.md` - these are snapshots, re-running reproduces
them exactly.)

## Generated SOPs

- `output/sops/SOP-001_HE_Staining_Kidney.md` - optimized staining
  incubation time and antibody/probe concentration, with rationale drawn
  directly from the screening data.
- `output/sops/SOP-002_Tissue_Sectioning_Kidney.md` - optimized
  sectioning thickness, with rationale (thickness drives both tissue
  integrity and multi-cycle marker-channel alignment, more severely for
  CODEX than Xenium).

## Project structure

```
simulate_protocol_experiments.py   # screening (60) + validation (60) experiments -> data/raw/
standardize.py                     # unified schema + spec flags -> data/processed/
analysis.py                        # Kruskal-Wallis, t-test, consistency metric -> output/
generate_sops.py                   # SOP generation from optimized params -> output/sops/
run_all.py                         # orchestrates the full pipeline
tests/test_pipeline.py             # pytest checks incl. best-vs-worst sanity check
```

## Tech stack

- Python 3
- pandas / numpy for data simulation and wrangling
- scipy.stats for Kruskal-Wallis and Welch's t-test
- pytest
