import os
import sys

import pandas as pd
import pytest

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, PROJECT_DIR)

import simulate_protocol_experiments as sim  # noqa: E402
import standardize  # noqa: E402
import analysis  # noqa: E402


@pytest.fixture(scope="module")
def pipeline_outputs():
    screening = sim.simulate_screening()
    best_variant_id, optimal_params = analysis.find_optimal_variant(screening)
    sim.simulate_validation(optimal_params)
    standardize.standardize()
    results = analysis.main()
    results["best_variant_id"] = best_variant_id
    results["optimal_params"] = optimal_params
    return results


def test_screening_has_expected_variants_and_count():
    screening = sim.simulate_screening()
    assert screening["variant_id"].nunique() == len(sim.SCREENING_VARIANTS)
    assert len(screening) == len(sim.SCREENING_VARIANTS) * len(sim.PLATFORMS) * sim.N_REPLICATES_SCREENING


def test_fifty_or_more_validation_experiments_on_optimized_protocol(pipeline_outputs):
    unified = pd.read_csv(os.path.join(standardize.PROCESSED_DIR, "unified_experiments.csv"))
    validation = unified[unified["stage"] == "validation"]
    assert len(validation) >= 50
    # all validation experiments use the single optimized protocol
    assert validation["staining_time_min"].nunique() == 1
    assert validation["antibody_conc_ug_ml"].nunique() == 1
    assert validation["thickness_um"].nunique() == 1


def test_total_experiments_at_least_fifty(pipeline_outputs):
    unified = pd.read_csv(os.path.join(standardize.PROCESSED_DIR, "unified_experiments.csv"))
    assert len(unified) >= 50


def test_image_quality_metrics_are_real_computed_columns(pipeline_outputs):
    unified = pd.read_csv(os.path.join(standardize.PROCESSED_DIR, "unified_experiments.csv"))
    for col in ["snr_db", "staining_uniformity_pct", "tissue_integrity_pct",
                "marker_alignment_pct", "image_quality_score"]:
        assert col in unified.columns
        assert unified[col].notna().all()
        assert unified[col].std() > 0  # not a constant/placeholder column


def test_kruskal_wallis_finds_significant_variant_differences(pipeline_outputs):
    assert pipeline_outputs["kruskal_p"] < 0.01


def test_optimized_variant_genuinely_outperforms_worst_variant(pipeline_outputs):
    """Injected sanity check: the chosen optimal protocol must be statistically
    and substantially better than the worst-performing screening variant, not
    just nominally ranked first."""
    assert pipeline_outputs["best_variant_id"] != pipeline_outputs["worst_variant_id"]
    ranked = pipeline_outputs["ranked"].set_index("variant_id")
    best_mean = ranked.loc[pipeline_outputs["best_variant_id"], "mean_quality"]
    worst_mean = ranked.loc[pipeline_outputs["worst_variant_id"], "mean_quality"]
    assert best_mean > worst_mean + 10  # substantial margin, not noise
    assert pipeline_outputs["ttest_p"] < 0.01


def test_consistency_metric_in_reasonable_band(pipeline_outputs):
    consistency = pipeline_outputs["overall_consistency"]
    assert 80.0 <= consistency <= 100.0


def test_sops_are_generated_with_real_content():
    import generate_sops
    he_path, sectioning_path = generate_sops.main()
    for path in (he_path, sectioning_path):
        assert os.path.exists(path)
        with open(path) as f:
            text = f.read()
        assert len(text) > 1000
        assert "Optimized" in text or "optimized" in text
        assert "ug/mL" in text or "um" in text
