"""
standardize.py

Concatenates the Stage 1 (screening) and Stage 2 (validation) experiment
tables from data/raw/ into one unified schema in data/processed/, and flags
every experiment against the fixed image-quality spec
(QUALITY_PASS_THRESHOLD from simulate_protocol_experiments.py).
"""

import os

import pandas as pd

from simulate_protocol_experiments import RAW_DIR, QUALITY_PASS_THRESHOLD

PROJECT_DIR = os.path.dirname(os.path.abspath(__file__))
PROCESSED_DIR = os.path.join(PROJECT_DIR, "data", "processed")


def standardize():
    screening_path = os.path.join(RAW_DIR, "screening_experiments.csv")
    validation_path = os.path.join(RAW_DIR, "validation_experiments.csv")
    frames = []
    for path in (screening_path, validation_path):
        if os.path.exists(path):
            frames.append(pd.read_csv(path))
    if not frames:
        raise FileNotFoundError("No raw experiment files found - run simulate_protocol_experiments.py first")

    unified = pd.concat(frames, ignore_index=True)
    unified["in_spec"] = unified["image_quality_score"] >= QUALITY_PASS_THRESHOLD

    os.makedirs(PROCESSED_DIR, exist_ok=True)
    out_path = os.path.join(PROCESSED_DIR, "unified_experiments.csv")
    unified.to_csv(out_path, index=False)

    print(f"Standardized {len(unified)} experiments -> {out_path}")
    print(f"  {unified['in_spec'].sum()} / {len(unified)} experiments in spec "
          f"(image_quality_score >= {QUALITY_PASS_THRESHOLD})")
    return unified


if __name__ == "__main__":
    standardize()
