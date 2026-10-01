# CODEX/Xenium Kidney Protocol Optimization - Analysis Report

Kruskal-Wallis across 10 screening variants: H = 49.026, p = 1.64e-07 (n=60)

Best variant: **thickness_3um** (mean image_quality_score = 91.35)
Worst variant: **time_60min** (mean image_quality_score = 53.28)
Best vs worst (Welch t-test): t = 41.337, p = 6.40e-11

Validation-stage image-quality consistency (optimized protocol, n=60): **91.67%** of experiments >= 87.5 quality-score threshold

## Per-platform consistency

| platform   |   n |   n_in_spec |   pct_consistent |   mean_quality |   std_quality |   cv_pct |
|:-----------|----:|------------:|-----------------:|---------------:|--------------:|---------:|
| CODEX      |  30 |          27 |            90    |          90.37 |          1.66 |     1.84 |
| Xenium     |  30 |          28 |            93.33 |          91.17 |          1.67 |     1.83 |
| overall    |  60 |          55 |            91.67 |          90.77 |          1.7  |     1.87 |
