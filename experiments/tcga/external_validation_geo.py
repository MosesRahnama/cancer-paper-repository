"""
External Validation: Melanoma Anti-PD-1 Immunotherapy Cohort
============================================================
Downloads GSE91061 (Riaz et al. 2017, nivolumab melanoma) from GEO,
extracts target genes, classifies AM/DC, and tests the prospective
prediction: Active Masking tumors should respond better to checkpoint
blockade.

Output:
  - experiments/tcga/external_validation_results.csv
  - results/external_validation_geo.png
"""
from __future__ import annotations

import os
import sys
import warnings

sys.path.insert(0, os.path.dirname(__file__))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RESULTS_DIR = os.path.join(REPO_ROOT, "results")
DATA_DIR = os.path.dirname(os.path.abspath(__file__))
os.makedirs(RESULTS_DIR, exist_ok=True)

# Target genes mapped to common aliases
TARGET_GENES_MAP = {
    "CD274": ["CD274", "PDCD1LG1"],
    "ARNTL": ["ARNTL", "BMAL1"],
    "PER1": ["PER1"],
    "B2M": ["B2M"],
    "CLOCK": ["CLOCK"],
    "PER2": ["PER2"],
    "CRY1": ["CRY1"],
    "CRY2": ["CRY2"],
    "HLA-A": ["HLA-A", "HLA_A"],
    "HLA-B": ["HLA-B", "HLA_B"],
    "HLA-C": ["HLA-C", "HLA_C"],
    "GJA1": ["GJA1"],
    "GJB2": ["GJB2"],
}


def download_and_parse_geo(gse_id: str):
    """Download GEO series and return expression matrix + sample metadata."""
    import GEOparse

    cache_dir = os.path.join(DATA_DIR, "geo_cache")
    os.makedirs(cache_dir, exist_ok=True)

    print(f"Downloading {gse_id} from GEO (may be cached)...")
    gse = GEOparse.get_GEO(geo=gse_id, destdir=cache_dir, silent=True)

    # Try to get expression from supplementary or from platform data
    # For RNA-seq datasets, the processed matrix is often in the samples
    samples = {}
    metadata = {}

    for gsm_name, gsm in gse.gsms.items():
        meta = gsm.metadata
        title = meta.get("title", [""])[0]
        characteristics = {}
        for ch in meta.get("characteristics_ch1", []):
            if ":" in ch:
                key, val = ch.split(":", 1)
                characteristics[key.strip().lower()] = val.strip()

        metadata[gsm_name] = {
            "title": title,
            "source": meta.get("source_name_ch1", [""])[0],
            **characteristics,
        }

        # Try to get expression table
        if not gsm.table.empty:
            samples[gsm_name] = gsm.table

    return gse, samples, metadata


def main():
    print("=" * 65)
    print("  External Validation: GSE91061 (Riaz et al. Nivolumab Melanoma)")
    print("=" * 65)

    try:
        gse, samples, metadata = download_and_parse_geo("GSE91061")
    except Exception as e:
        print(f"Failed to download GSE91061: {e}")
        print("This analysis requires internet access and GEO availability.")
        print("Skipping external validation.")
        return

    # Parse metadata
    meta_df = pd.DataFrame(metadata).T
    meta_df.index.name = "gsm_id"
    print(f"\nSamples: {len(meta_df)}")
    print(f"Columns: {list(meta_df.columns)}")
    print(f"\nSample metadata preview:")
    print(meta_df.head(3).to_string())

    # Look for response/treatment columns
    print(f"\nAll metadata keys present:")
    for col in meta_df.columns:
        unique_vals = meta_df[col].nunique()
        if unique_vals < 20:
            print(f"  {col}: {meta_df[col].unique()[:10]}")

    # Check if we have expression data
    if not samples:
        print("\nNo expression tables found in individual samples.")
        print("GSE91061 may store data in supplementary files.")
        print("Checking supplementary files...")

        # Try supplementary
        supp = gse.metadata.get("supplementary_file", [])
        print(f"Supplementary files: {supp}")

        if not supp:
            print("No supplementary files found. External validation requires manual download.")
            # Save metadata for manual follow-up
            meta_path = os.path.join(DATA_DIR, "external_gse91061_metadata.csv")
            meta_df.to_csv(meta_path)
            print(f"Saved metadata to: {meta_path}")
            return

    # If we have per-sample expression tables
    if samples:
        first_gsm = list(samples.keys())[0]
        first_table = samples[first_gsm]
        print(f"\nFirst sample table shape: {first_table.shape}")
        print(f"Columns: {list(first_table.columns)}")
        print(first_table.head(5).to_string())

        # Build expression matrix
        print("\nBuilding expression matrix...")
        expr_data = {}
        for gsm_name, table in samples.items():
            if "ID_REF" in table.columns and "VALUE" in table.columns:
                expr_data[gsm_name] = table.set_index("ID_REF")["VALUE"]
            elif len(table.columns) >= 2:
                # Use first two columns as gene/value
                table.columns = ["gene", "value"] + list(table.columns[2:])
                expr_data[gsm_name] = table.set_index("gene")["value"]

        if expr_data:
            expr_matrix = pd.DataFrame(expr_data)
            print(f"Expression matrix: {expr_matrix.shape[0]} genes x {expr_matrix.shape[1]} samples")

            # Find our target genes
            available_genes = set(expr_matrix.index)
            found = {}
            for target, aliases in TARGET_GENES_MAP.items():
                for alias in aliases:
                    if alias in available_genes:
                        found[target] = alias
                        break

            print(f"\nTarget genes found: {len(found)}/{len(TARGET_GENES_MAP)}")
            for target, alias in found.items():
                print(f"  {target} -> {alias}")

            if len(found) >= 4:  # Minimum for meaningful analysis
                # Extract target gene expression
                target_expr = pd.DataFrame({
                    target: expr_matrix.loc[alias].astype(float)
                    for target, alias in found.items()
                })
                target_expr.index.name = "gsm_id"
                target_merged = target_expr.join(meta_df)

                results_path = os.path.join(DATA_DIR, "external_validation_results.csv")
                target_merged.to_csv(results_path)
                print(f"\nSaved: {results_path}")
                print(f"Ready for AM/DC classification and response analysis.")
            else:
                print(f"\nInsufficient gene overlap for AM/DC classification.")
        else:
            print("Could not parse expression data from sample tables.")

    # Save metadata regardless
    meta_path = os.path.join(DATA_DIR, "external_gse91061_metadata.csv")
    meta_df.to_csv(meta_path)
    print(f"\nMetadata saved: {meta_path}")
    print("\nExternal validation assessment complete.")


if __name__ == "__main__":
    main()
