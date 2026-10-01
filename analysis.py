"""
analysis.py

Real statistical comparison of the Stage 1 screening protocol variants,
selection of the optimal protocol, and a real image-quality-consistency
metric computed from the Stage 2 validation experiments run with that
optimized protocol.

Outputs (written under output/):
  variant_comparison.csv   per-variant mean/std image_quality_score + params
  consistency_summary.csv  validation-stage consistency metric, by platform and overall
  analysis_report.md       human-readable summary of the whole analysis
"""

import os

import numpy as np
import pandas as pd
from scipy import stats

from simulate_protocol_experiments import QUALITY_PASS_THRESHOLD

PROJECT_DIR = os.path.dirname(os.path.abspath(__file__))
PROCESSED_DIR = os.path.join(PROJECT_DIR, "data", "processed")
OUTPUT_DIR = os.path.join(PROJECT_DIR, "output")


def load_unified():
    return pd.read_csv(os.path.join(PROCESSED_DIR, "unified_experiments.csv"))


def variant_stats(screening):
    """Mean/std image_quality_score per screening variant, plus its protocol params."""
    grouped = screening.groupby("variant_id").agg(
        n=("image_quality_score", "size"),
        mean_quality=("image_quality_score", "mean"),
        std_quality=("image_quality_score", "std"),
        staining_time_min=("staining_time_min", "first"),
        antibody_conc_ug_ml=("antibody_conc_ug_ml", "first"),
        thickness_um=("thickness_um", "first"),
    ).reset_index().sort_values("mean_quality", ascending=False)
    return grouped


def kruskal_across_variants(screening):
    groups = [g["image_quality_score"].values for _, g in screening.groupby("variant_id")]
    h, p = stats.kruskal(*groups)
    return h, p


def best_vs_worst_ttest(screening, best_variant_id, worst_variant_id):
    best = screening.loc[screening["variant_id"] == best_variant_id, "image_quality_score"]
    worst = screening.loc[screening["variant_id"] == worst_variant_id, "image_quality_score"]
    t, p = stats.ttest_ind(best, worst, equal_var=False)
    return t, p, best.mean(), worst.mean()


def find_optimal_variant(screening):
    """Returns (variant_id, params dict) for the highest-mean-quality screening variant."""
    ranked = variant_stats(screening)
    top = ranked.iloc[0]
    params = dict(
        staining_time_min=int(top["staining_time_min"]),
        antibody_conc_ug_ml=float(top["antibody_conc_ug_ml"]),
        thickness_um=float(top["thickness_um"]),
    )
    return top["variant_id"], params


def consistency_metric(validation):
    """Image-quality-consistency: % of validation experiments meeting the quality spec.

    Computed overall and per platform. Also reports the coefficient of
    variation of image_quality_score as a secondary consistency measure
    (lower CV = more consistent).
    """
    rows = []
    for platform, g in list(validation.groupby("platform")) + [("overall", validation)]:
        n = len(g)
        n_pass = int((g["image_quality_score"] >= QUALITY_PASS_THRESHOLD).sum())
        pct_consistent = 100.0 * n_pass / n
        cv_pct = 100.0 * g["image_quality_score"].std() / g["image_quality_score"].mean()
        rows.append(dict(
            platform=platform, n=n, n_in_spec=n_pass,
            pct_consistent=round(pct_consistent, 2),
            mean_quality=round(g["image_quality_score"].mean(), 2),
            std_quality=round(g["image_quality_score"].std(), 2),
            cv_pct=round(cv_pct, 2),
        ))
    return pd.DataFrame(rows)


def main():
    unified = load_unified()
    screening = unified[unified["stage"] == "screening"].copy()
    validation = unified[unified["stage"] == "validation"].copy()

    ranked = variant_stats(screening)
    h, p_kruskal = kruskal_across_variants(screening)
    best_variant_id = ranked.iloc[0]["variant_id"]
    worst_variant_id = ranked.iloc[-1]["variant_id"]
    t, p_ttest, best_mean, worst_mean = best_vs_worst_ttest(screening, best_variant_id, worst_variant_id)

    consistency = consistency_metric(validation)
    overall_consistency = consistency.loc[consistency["platform"] == "overall", "pct_consistent"].iloc[0]

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    ranked.to_csv(os.path.join(OUTPUT_DIR, "variant_comparison.csv"), index=False)
    consistency.to_csv(os.path.join(OUTPUT_DIR, "consistency_summary.csv"), index=False)

    report_lines = [
        "# CODEX/Xenium Kidney Protocol Optimization - Analysis Report",
        "",
        f"Kruskal-Wallis across {ranked.shape[0]} screening variants: "
        f"H = {h:.3f}, p = {p_kruskal:.2e} (n={len(screening)})",
        "",
        f"Best variant: **{best_variant_id}** (mean image_quality_score = {best_mean:.2f})",
        f"Worst variant: **{worst_variant_id}** (mean image_quality_score = {worst_mean:.2f})",
        f"Best vs worst (Welch t-test): t = {t:.3f}, p = {p_ttest:.2e}",
        "",
        f"Validation-stage image-quality consistency (optimized protocol, "
        f"n={len(validation)}): **{overall_consistency:.2f}%** of experiments >= "
        f"{QUALITY_PASS_THRESHOLD} quality-score threshold",
        "",
        "## Per-platform consistency",
        "",
        consistency.to_markdown(index=False),
    ]
    with open(os.path.join(OUTPUT_DIR, "analysis_report.md"), "w") as f:
        f.write("\n".join(report_lines) + "\n")

    print(f"Kruskal-Wallis across variants: H = {h:.3f}, p = {p_kruskal:.2e}")
    print(f"Best variant: {best_variant_id} (mean = {best_mean:.2f}) vs "
          f"worst: {worst_variant_id} (mean = {worst_mean:.2f}); "
          f"t = {t:.3f}, p = {p_ttest:.2e}")
    print(f"Validation-stage image-quality consistency: {overall_consistency:.2f}% "
          f"(n={len(validation)})")

    return dict(
        ranked=ranked, kruskal_h=h, kruskal_p=p_kruskal,
        best_variant_id=best_variant_id, worst_variant_id=worst_variant_id,
        ttest_t=t, ttest_p=p_ttest, consistency=consistency,
        overall_consistency=overall_consistency,
    )


if __name__ == "__main__":
    main()
