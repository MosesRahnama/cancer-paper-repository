"""
TCGA Expression Data Extraction Pipeline
=========================================
Queries ISB-CGC BigQuery for TCGA-SKCM and TCGA-LUAD expression data,
extracts target genes for Hypothesis 1 validation, and saves CSV outputs.
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))

import pandas as pd

from tcga_config import (
    DATA_DIR,
    GCS_PREFIX,
    get_bq_client,
    get_gcs_bucket,
    safe_alias,
)

TARGET_GENES = [
    "CD274",
    "PDCD1LG2",
    "PDCD1",
    "HLA-A",
    "HLA-B",
    "HLA-C",
    "B2M",
    "GJA1",
    "GJB2",
    "GJA5",
    "GJB6",
    "ARNTL",
    "CLOCK",
    "PER1",
    "PER2",
    "CRY1",
    "CRY2",
    "MYC",
    "TP53",
    "CDH1",
    "VIM",
]

TARGET_PROJECTS = ["TCGA-SKCM", "TCGA-LUAD"]


def upload_df(bucket, gcs_path: str, df: pd.DataFrame):
    blob = bucket.blob(gcs_path)
    blob.upload_from_string(df.to_csv(index=False), content_type="text/csv")
    print(f"  Uploaded: gs://{bucket.name}/{gcs_path} ({len(df)} rows)")


def main():
    parser = argparse.ArgumentParser(description="Extract TCGA expression matrices from BigQuery.")
    parser.add_argument(
        "--output-dir",
        default=DATA_DIR,
        help="Local output directory for extracted CSV files.",
    )
    parser.add_argument(
        "--no-gcs",
        action="store_true",
        help="Skip upload to GCS and only write local CSV files.",
    )
    args = parser.parse_args()

    bq_client = get_bq_client()
    bucket = None
    if not args.no_gcs:
        try:
            bucket = get_gcs_bucket()
            if not bucket.exists():
                raise RuntimeError(f"GCS bucket {bucket.name} does not exist.")
            print(f"Verified GCS bucket gs://{bucket.name}")
        except Exception as e:
            print(f"GCS setup failed; continuing with local outputs only: {e}")
            bucket = None

    gene_list_sql = ", ".join(f"'{g}'" for g in TARGET_GENES)
    project_list_sql = ", ".join(f"'{p}'" for p in TARGET_PROJECTS)

    query = f"""
    SELECT
      project_short_name,
      case_barcode,
      sample_barcode,
      sample_type_name,
      {', '.join(
        f"MAX(IF(gene_name = '{g}', tpm_unstranded, NULL)) AS `{safe_alias(g)}`"
        for g in TARGET_GENES
      )}
    FROM `isb-cgc-bq.TCGA.RNAseq_hg38_gdc_current`
    WHERE project_short_name IN ({project_list_sql})
      AND gene_name IN ({gene_list_sql})
    GROUP BY project_short_name, case_barcode, sample_barcode, sample_type_name
    ORDER BY project_short_name, case_barcode
    """

    print(f"\nQuerying BigQuery for {len(TARGET_GENES)} genes across {TARGET_PROJECTS}...")
    df = bq_client.query(query).to_dataframe()
    print(f"  Retrieved {len(df)} sample rows.")
    print(f"  Projects: {df['project_short_name'].value_counts().to_dict()}")
    print(f"  Sample types: {df['sample_type_name'].value_counts().to_dict()}")

    query_counts = f"""
    SELECT
      project_short_name,
      case_barcode,
      sample_barcode,
      sample_type_name,
      {', '.join(
        f"MAX(IF(gene_name = '{g}', unstranded, NULL)) AS `{safe_alias(g)}_counts`"
        for g in TARGET_GENES
      )}
    FROM `isb-cgc-bq.TCGA.RNAseq_hg38_gdc_current`
    WHERE project_short_name IN ({project_list_sql})
      AND gene_name IN ({gene_list_sql})
    GROUP BY project_short_name, case_barcode, sample_barcode, sample_type_name
    ORDER BY project_short_name, case_barcode
    """
    print("\nQuerying raw counts...")
    df_counts = bq_client.query(query_counts).to_dataframe()
    print(f"  Retrieved {len(df_counts)} sample rows (counts).")

    os.makedirs(args.output_dir, exist_ok=True)
    out_tpm = os.path.join(args.output_dir, "tcga_expression_tpm.csv")
    out_counts = os.path.join(args.output_dir, "tcga_expression_counts.csv")
    df.to_csv(out_tpm, index=False)
    df_counts.to_csv(out_counts, index=False)
    print(f"\nSaved local outputs:\n  {out_tpm}\n  {out_counts}")

    if bucket is not None:
        print(f"\nUploading to gs://{bucket.name}/{GCS_PREFIX}/data/...")
        upload_df(bucket, f"{GCS_PREFIX}/data/tcga_expression_tpm.csv", df)
        upload_df(bucket, f"{GCS_PREFIX}/data/tcga_expression_counts.csv", df_counts)

    print("\n" + "=" * 60)
    print("EXTRACTION COMPLETE")
    print("=" * 60)
    for proj in TARGET_PROJECTS:
        subset = df[df["project_short_name"] == proj]
        tumor = subset[subset["sample_type_name"].str.contains("Tumor", case=False, na=False)]
        normal = subset[~subset["sample_type_name"].str.contains("Tumor", case=False, na=False)]
        print(f"\n  {proj}:")
        print(f"    Total samples: {len(subset)}")
        print(f"    Tumor samples: {len(tumor)}")
        print(f"    Normal/other:  {len(normal)}")
        safe_cols = [safe_alias(g) for g in TARGET_GENES]
        available_cols = [c for c in safe_cols if c in subset.columns]
        gene_nulls = subset[available_cols].isnull().sum()
        genes_with_data = (gene_nulls < len(subset)).sum()
        print(f"    Genes with data: {genes_with_data}/{len(TARGET_GENES)}")

    if bucket is not None:
        print(f"\n  GCS location: gs://{bucket.name}/{GCS_PREFIX}/data/")
    print("  Ready for correlation analysis.\n")


if __name__ == "__main__":
    main()
