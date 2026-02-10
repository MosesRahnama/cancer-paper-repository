"""
Convert Thorsson et al. 2018 Immunity supplemental data (mmc2.xlsx)
to the purity/immune covariate CSV format for robustness_check.py
"""
import os
import pandas as pd

# Source file from Downloads
SOURCE_FILE = r"C:\Users\Moses\OpComp\Mina Analytics Dropbox\Moses Rahnama\PC\Downloads\mmc2.xlsx"

# Target location in cancer-paper-repository
TARGET_FILE = r"C:\Users\Moses\Cancer\cancer-paper-repository\experiments\tcga\tcga_purity_immune_covariates.csv"

# Cancer types in our analysis (abbreviated format used in Thorsson data)
TARGET_PROJECTS = [
    "SKCM",  # Melanoma
    "LUAD",  # Lung adenocarcinoma
    "BRCA",  # Breast invasive carcinoma
    "COAD",  # Colon adenocarcinoma
    "HNSC",  # Head and neck squamous cell carcinoma
    "LUSC",  # Lung squamous cell carcinoma
]

def main():
    print("Loading Thorsson Immunity data from mmc2.xlsx...")

    # Read the Excel file (should be first sheet with sample-level data)
    df = pd.read_excel(SOURCE_FILE, sheet_name=0)

    print(f"  Loaded {len(df)} rows")
    print(f"\n  Available columns:")
    for col in df.columns:
        print(f"    - {col}")

    # Identify the key columns
    # Expected: TCGA Participant Barcode, TCGA Study, Leukocyte Fraction,
    #           Lymphocyte Infiltration Signature Score, and potentially CPE/purity columns

    # Map to our expected column names
    column_mapping = {}

    # Case barcode (patient ID)
    for variant in ["TCGA Participant Barcode", "TCGA Patient Barcode", "Participant Barcode", "barcode"]:
        if variant in df.columns:
            column_mapping[variant] = "case_barcode"
            break

    # Study/project
    for variant in ["TCGA Study", "Study", "project_id"]:
        if variant in df.columns:
            column_mapping[variant] = "project_short_name"
            break

    # Leukocyte/immune infiltration
    for variant in ["Leukocyte Fraction", "Leukocyte Score"]:
        if variant in df.columns:
            column_mapping[variant] = "leukocyte_fraction"
            break

    # Lymphocyte infiltration
    for variant in ["Lymphocyte Infiltration Signature Score", "Lymphocyte Score"]:
        if variant in df.columns:
            column_mapping[variant] = "lymphocyte_infiltration"
            break

    # Stromal fraction (related to tumor microenvironment)
    for variant in ["Stromal Fraction"]:
        if variant in df.columns:
            column_mapping[variant] = "stromal_fraction"
            break

    # IFN-gamma response (immune activation marker)
    for variant in ["IFN-gamma Response"]:
        if variant in df.columns:
            column_mapping[variant] = "ifn_gamma_response"
            break

    # TGF-beta response (immunosuppression marker)
    for variant in ["TGF-beta Response"]:
        if variant in df.columns:
            column_mapping[variant] = "tgfb_response"
            break

    # Macrophage regulation
    for variant in ["Macrophage Regulation"]:
        if variant in df.columns:
            column_mapping[variant] = "macrophage_regulation"
            break

    # Immune subtype
    for variant in ["Immune Subtype"]:
        if variant in df.columns:
            column_mapping[variant] = "immune_subtype"
            break

    # Tumor purity (check for CPE, ESTIMATE, etc.)
    for variant in ["CPE", "ESTIMATE Score", "Purity", "Tumor Purity"]:
        if variant in df.columns:
            column_mapping[variant] = "purity"
            break

    print(f"\n  Column mapping:")
    for old, new in column_mapping.items():
        print(f"    {old} -> {new}")

    if "case_barcode" not in column_mapping.values():
        raise ValueError("Could not find patient barcode column!")

    # Rename columns
    df_out = df.rename(columns=column_mapping)

    # Keep only mapped columns
    keep_cols = [col for col in df_out.columns if col in column_mapping.values()]
    df_out = df_out[keep_cols]

    # Check what project values are actually in the data
    if "project_short_name" in df_out.columns:
        print(f"\n  Unique project values in data:")
        unique_projects = df_out["project_short_name"].unique()
        for proj in sorted(unique_projects)[:20]:  # Show first 20
            count = (df_out["project_short_name"] == proj).sum()
            print(f"    {proj}: {count}")
        if len(unique_projects) > 20:
            print(f"    ... and {len(unique_projects) - 20} more")

        # Filter to our 6 cancer types
        print(f"\n  Filtering to {len(TARGET_PROJECTS)} cancer types...")
        print(f"  Before filter: {len(df_out)} rows")
        df_out = df_out[df_out["project_short_name"].isin(TARGET_PROJECTS)]
        print(f"  After filter: {len(df_out)} rows")

        if len(df_out) > 0:
            print(f"\n  Breakdown by cancer type:")
            for proj in TARGET_PROJECTS:
                count = (df_out["project_short_name"] == proj).sum()
                if count > 0:
                    print(f"    {proj}: {count}")

    # Clean case_barcode (remove -01A, -11A suffixes to match clinical data format)
    if "case_barcode" in df_out.columns:
        df_out["case_barcode"] = df_out["case_barcode"].str[:12]  # Keep only TCGA-XX-XXXX format

    # Add TCGA- prefix to project names to match format used in main analysis
    if "project_short_name" in df_out.columns:
        df_out["project_short_name"] = "TCGA-" + df_out["project_short_name"]

    # Save to target location
    os.makedirs(os.path.dirname(TARGET_FILE), exist_ok=True)
    df_out.to_csv(TARGET_FILE, index=False)

    print(f"\nSaved to: {TARGET_FILE}")
    print(f"  {len(df_out)} samples across {len(TARGET_PROJECTS)} cancer types")
    print(f"  Columns: {', '.join(df_out.columns)}")

    # Show summary statistics
    print(f"\n  Summary statistics:")
    print(df_out.describe())


if __name__ == "__main__":
    main()
