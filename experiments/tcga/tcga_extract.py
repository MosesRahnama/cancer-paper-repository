"""
TCGA Expression Data Extraction Pipeline
=========================================
Queries ISB-CGC BigQuery for TCGA-SKCM (Melanoma) and TCGA-LUAD (Lung Adenocarcinoma)
gene expression data, extracts target genes for Hypothesis 1 validation
(PD-L1 vs gap junctions vs circadian coherence), and exports results to GCS.

Target genes:
  - Immune checkpoint: CD274 (PD-L1), PDCD1LG2 (PD-L2), PDCD1 (PD-1)
  - MHC-I antigen presentation: HLA-A, HLA-B, HLA-C, B2M
  - Gap junctions: GJA1 (Cx43), GJB2 (Cx26), GJA5 (Cx40), GJB6 (Cx30)
  - Circadian clock: ARNTL (BMAL1), CLOCK, PER1, PER2, CRY1, CRY2
  - Stemness/entropy: MYC, TP53, CDH1 (E-cadherin), VIM (Vimentin)
"""
import os
import sys

os.environ["GOOGLE_APPLICATION_CREDENTIALS"] = (
    r"C:\Users\Moses\commercials\infra\secrets\gcp-credentials.json"
)

from google.cloud import bigquery, storage
import pandas as pd
import io

# ── Configuration ──────────────────────────────────────────────────────────
PROJECT_ID = "kernel-o6"
GCS_BUCKET = "kernel-o6"
GCS_PREFIX = "tcga_cancer_boundary/data"

TARGET_GENES = [
    # Immune checkpoint
    "CD274", "PDCD1LG2", "PDCD1",
    # MHC-I
    "HLA-A", "HLA-B", "HLA-C", "B2M",
    # Gap junctions (connexins)
    "GJA1", "GJB2", "GJA5", "GJB6",
    # Circadian clock
    "ARNTL", "CLOCK", "PER1", "PER2", "CRY1", "CRY2",
    # Differentiation / stemness markers
    "MYC", "TP53", "CDH1", "VIM",
]

TARGET_PROJECTS = ["TCGA-SKCM", "TCGA-LUAD"]


def safe_alias(gene_name):
    """Convert gene name to valid SQL alias (replace hyphens with underscores)."""
    return gene_name.replace("-", "_")


def main():
    bq_client = bigquery.Client(project=PROJECT_ID)
    gcs_client = storage.Client(project=PROJECT_ID)
    bucket = gcs_client.bucket(GCS_BUCKET)

    # ── Step 1: Check bucket access ────────────────────────────────────────
    print(f"Verifying GCS bucket gs://{GCS_BUCKET}...")
    if bucket.exists():
        print(f"  Bucket exists.")
    else:
        print(f"  Bucket not found. Creating gs://{GCS_BUCKET}...")
        bucket = gcs_client.create_bucket(GCS_BUCKET, location="US")
        print(f"  Created bucket gs://{GCS_BUCKET}")

    # ── Step 2: Query expression data (pivoted wide) ───────────────────────
    gene_list_sql = ", ".join(f"'{g}'" for g in TARGET_GENES)
    project_list_sql = ", ".join(f"'{p}'" for p in TARGET_PROJECTS)

    # Query: one row per sample, one column per gene (TPM values)
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
    print(f"  This queries ISB-CGC public dataset (no download needed).")

    df = bq_client.query(query).to_dataframe()

    print(f"  Retrieved {len(df)} sample rows.")
    print(f"  Projects: {df['project_short_name'].value_counts().to_dict()}")
    print(f"  Sample types: {df['sample_type_name'].value_counts().to_dict()}")
    print(f"  Columns: {list(df.columns)}")
    print(f"\n  First 3 rows:")
    print(df.head(3).to_string(index=False))

    # ── Step 3: Also get raw counts for robust analysis ────────────────────
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

    print(f"\nQuerying raw counts...")
    df_counts = bq_client.query(query_counts).to_dataframe()
    print(f"  Retrieved {len(df_counts)} sample rows (counts).")

    # ── Step 4: Upload to GCS ──────────────────────────────────────────────
    print(f"\nUploading to gs://{GCS_BUCKET}/{GCS_PREFIX}/...")

    # TPM expression matrix
    csv_tpm = df.to_csv(index=False)
    blob_tpm = bucket.blob(f"{GCS_PREFIX}/tcga_expression_tpm.csv")
    blob_tpm.upload_from_string(csv_tpm, content_type="text/csv")
    print(f"  Uploaded: gs://{GCS_BUCKET}/{GCS_PREFIX}/tcga_expression_tpm.csv ({len(df)} rows)")

    # Raw counts matrix
    csv_counts = df_counts.to_csv(index=False)
    blob_counts = bucket.blob(f"{GCS_PREFIX}/tcga_expression_counts.csv")
    blob_counts.upload_from_string(csv_counts, content_type="text/csv")
    print(f"  Uploaded: gs://{GCS_BUCKET}/{GCS_PREFIX}/tcga_expression_counts.csv ({len(df_counts)} rows)")

    # Also save locally for immediate analysis
    local_dir = r"C:\Users\Moses\Cancer\data"
    os.makedirs(local_dir, exist_ok=True)
    df.to_csv(os.path.join(local_dir, "tcga_expression_tpm.csv"), index=False)
    df_counts.to_csv(os.path.join(local_dir, "tcga_expression_counts.csv"), index=False)
    print(f"  Also saved locally to {local_dir}")

    # ── Step 5: Summary statistics ─────────────────────────────────────────
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

        # Check gene coverage
        safe_cols = [safe_alias(g) for g in TARGET_GENES]
        available_cols = [c for c in safe_cols if c in subset.columns]
        gene_nulls = subset[available_cols].isnull().sum()
        genes_with_data = (gene_nulls < len(subset)).sum()
        print(f"    Genes with data: {genes_with_data}/{len(TARGET_GENES)}")

    print(f"\n  GCS location: gs://{GCS_BUCKET}/{GCS_PREFIX}/")
    print(f"  Ready for correlation analysis.\n")


if __name__ == "__main__":
    main()
