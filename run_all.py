"""
run_all.py

Runs the full pipeline end-to-end:
  1. Simulate Stage 1 screening experiments (protocol-variant sweep) -> data/raw/
  2. Identify the optimal protocol variant from the screening data
  3. Simulate Stage 2 validation experiments using the optimized protocol -> data/raw/
  4. Standardize screening + validation into one unified schema -> data/processed/
  5. Run statistical comparisons and compute the image-quality-consistency metric -> output/
  6. Generate the H&E staining and tissue sectioning SOPs -> output/sops/

Run with:  python run_all.py
"""

import simulate_protocol_experiments as sim
import standardize
import analysis
import generate_sops


def main():
    print("--- Stage 1: simulating screening experiments ---")
    screening = sim.simulate_screening()

    print("\n--- Identifying optimal protocol variant ---")
    best_variant_id, optimal_params = analysis.find_optimal_variant(screening)
    print(f"Optimal variant (provisional, from screening means): {best_variant_id} -> {optimal_params}")

    print("\n--- Stage 2: simulating validation experiments with optimized protocol ---")
    sim.simulate_validation(optimal_params)

    print("\n--- Standardizing into unified schema ---")
    standardize.standardize()

    print("\n--- Statistical analysis + consistency metric ---")
    analysis.main()

    print("\n--- Generating SOPs ---")
    generate_sops.main()


if __name__ == "__main__":
    main()
