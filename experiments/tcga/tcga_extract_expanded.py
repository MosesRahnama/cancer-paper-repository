"""
Expanded TCGA extraction: 6 cancer types + clinical + diagnosis data.
Queries ISB-CGC BigQuery and exports to GCS + local CSV.
"""
import argparse
import sys
import os
sys.path.insert(0, os.path.dirname(__file__))

from tcga_config import (
    get_bq_client, get_gcs_bucket, safe_alias,
    TARGET_GENES_RAW, ALL_PROJECTS, DATA_DIR, GCS_PREFIX,
    TCGA_DEFAULT_RELEASE_TAG, resolve_tcga_release_tag, tcga_table,
    write_extraction_manifest,
)
import pandas as pd


def extract_expression(client, expression_table):
    """Extract pivoted expression data for all 6 cancer types."""
    gene_list = ", ".join(f"'{g}'" for g in TARGET_GENES_RAW)
    proj_list = ", ".join(f"'{p}'" for p in ALL_PROJECTS)

    pivot_cols = ", ".join(
        f"MAX(IF(gene_name = '{g}', tpm_unstranded, NULL)) AS `{safe_alias(g)}`"
        for g in TARGET_GENES_RAW
    )

    query = f"""
    SELECT
      project_short_name,
      case_barcode,
      sample_barcode,
      sample_type_name,
      {pivot_cols}
    FROM `{expression_table}`
    WHERE project_short_name IN ({proj_list})
      AND gene_name IN ({gene_list})
    GROUP BY project_short_name, case_barcode, sample_barcode, sample_type_name
    ORDER BY project_short_name, case_barcode
    """

    print("Querying expression data for 6 cancer types...")
    df = client.query(query).to_dataframe()
    print(f"  Retrieved {len(df)} samples.")
    print(f"  Breakdown: {df['project_short_name'].value_counts().to_dict()}")
    return df


def extract_clinical(client, clinical_table, diagnosis_table):
    """Extract clinical + diagnosis data joined on case_id."""
    proj_list = ", ".join(f"'{p}'" for p in ALL_PROJECTS)

    query = f"""
    SELECT
      c.submitter_id AS case_barcode,
      c.proj__project_id AS project_short_name,
      c.primary_site,
      c.demo__vital_status AS vital_status,
      c.demo__days_to_death AS days_to_death,
      c.demo__gender AS gender,
      c.demo__age_at_index AS age_at_index,
      d.diag__ajcc_pathologic_stage AS ajcc_stage,
      d.diag__tumor_grade AS tumor_grade,
      d.diag__days_to_last_follow_up AS days_to_last_follow_up,
      d.diag__age_at_diagnosis AS age_at_diagnosis,
      d.diag__primary_diagnosis AS primary_diagnosis
    FROM `{clinical_table}` c
    LEFT JOIN `{diagnosis_table}` d
      ON c.case_id = d.case_id
    WHERE c.proj__project_id IN ({proj_list})
    """

    print("\nQuerying clinical + diagnosis data...")
    df = client.query(query).to_dataframe()

    # Deduplicate: some patients have multiple diagnosis records
    # Keep the first (primary) diagnosis per patient
    df = df.sort_values("case_barcode").drop_duplicates(subset=["case_barcode"], keep="first")

    print(f"  Retrieved {len(df)} clinical records (deduplicated).")
    if "vital_status" in df.columns:
        print(f"  Vital status: {df['vital_status'].value_counts().to_dict()}")
    if "ajcc_stage" in df.columns:
        non_null = df["ajcc_stage"].notna().sum()
        print(f"  AJCC stage available: {non_null}/{len(df)}")
    return df


def main():
    parser = argparse.ArgumentParser(
        description="Expanded TCGA extraction with deterministic table-version locking."
    )
    parser.add_argument(
        "--output-dir",
        default=DATA_DIR,
        help="Directory for output CSV and manifest files.",
    )
    parser.add_argument(
        "--release-tag",
        default=TCGA_DEFAULT_RELEASE_TAG,
        help="TCGA table release tag (default pinned): rNN or current.",
    )
    parser.add_argument(
        "--allow-current",
        action="store_true",
        help="Allow non-deterministic *_current table extraction (explicit opt-in).",
    )
    parser.add_argument(
        "--no-gcs",
        action="store_true",
        help="Skip GCS upload and write local files only.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Write extraction manifest only; do not query BigQuery.",
    )
    args = parser.parse_args()

    release_tag = resolve_tcga_release_tag(
        release_tag=args.release_tag,
        allow_current=args.allow_current,
    )
    tables = {
        "expression": tcga_table("RNAseq_hg38_gdc", release_tag=release_tag, allow_current=args.allow_current),
        "clinical": tcga_table("clinical_gdc", release_tag=release_tag, allow_current=args.allow_current),
        "diagnosis": tcga_table("clinical_diagnosis_gdc", release_tag=release_tag, allow_current=args.allow_current),
    }

    os.makedirs(args.output_dir, exist_ok=True)
    manifest_path = os.path.join(args.output_dir, "tcga_extract_expanded_manifest.json")
    write_extraction_manifest(
        output_path=manifest_path,
        script_name="tcga_extract_expanded.py",
        release_tag=release_tag,
        tables=tables,
        projects=ALL_PROJECTS,
        target_genes=TARGET_GENES_RAW,
    )

    print(f"Resolved TCGA release tag: {release_tag}")
    print(f"Expression table: {tables['expression']}")
    print(f"Clinical table:   {tables['clinical']}")
    print(f"Diagnosis table:  {tables['diagnosis']}")
    print(f"Wrote manifest:   {manifest_path}")

    if args.dry_run:
        print("Dry-run complete. No BigQuery queries were executed.")
        return

    client = get_bq_client()
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

    # ── Expression data ────────────────────────────────────────────────────
    expr_path = os.path.join(args.output_dir, "tcga_expanded_tpm.csv")
    if os.path.exists(expr_path):
        print(f"Expression data already exists at {expr_path}, skipping re-extraction.")
        df_expr = pd.read_csv(expr_path)
    else:
        df_expr = extract_expression(client, tables["expression"])
        df_expr.to_csv(expr_path, index=False)
        print(f"\n  Saved: {expr_path}")
        if bucket is not None:
            blob = bucket.blob(f"{GCS_PREFIX}/data/tcga_expanded_tpm.csv")
            blob.upload_from_filename(expr_path)
            print(f"  -> gs://{bucket.name}/{GCS_PREFIX}/data/tcga_expanded_tpm.csv")

    # ── Clinical data ──────────────────────────────────────────────────────
    df_clin = extract_clinical(client, tables["clinical"], tables["diagnosis"])
    clin_path = os.path.join(args.output_dir, "tcga_clinical.csv")
    df_clin.to_csv(clin_path, index=False)
    print(f"\n  Saved: {clin_path}")

    if bucket is not None:
        blob = bucket.blob(f"{GCS_PREFIX}/data/tcga_clinical.csv")
        blob.upload_from_filename(clin_path)
        print(f"  -> gs://{bucket.name}/{GCS_PREFIX}/data/tcga_clinical.csv")

    # ── Summary ────────────────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("EXPANDED EXTRACTION COMPLETE")
    print("=" * 60)
    for proj in ALL_PROJECTS:
        sub = df_expr[df_expr["project_short_name"] == proj]
        tumor = sub[sub["sample_type_name"].str.contains("Tumor|Metastatic", case=False, na=False)]
        normal = sub[~sub["sample_type_name"].str.contains("Tumor|Metastatic", case=False, na=False)]
        clin_sub = df_clin[df_clin["project_short_name"] == proj]
        alive = (clin_sub["vital_status"] == "Alive").sum() if "vital_status" in clin_sub.columns else "?"
        dead = (clin_sub["vital_status"] == "Dead").sum() if "vital_status" in clin_sub.columns else "?"
        stage_avail = clin_sub["ajcc_stage"].notna().sum() if "ajcc_stage" in clin_sub.columns else "?"
        print(f"  {proj:12s}: {len(tumor):4d} tumor, {len(normal):3d} normal | {alive} alive, {dead} dead | {stage_avail} staged")

    print(f"\n  Total expression samples: {len(df_expr)}")
    print(f"  Total clinical records: {len(df_clin)}")


if __name__ == "__main__":
    main()
