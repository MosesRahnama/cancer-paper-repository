"""
Shared configuration and utilities for TCGA boundary-logic analysis.
"""
import os
import numpy as np
import pandas as pd

# ── Paths ──────────────────────────────────────────────────────────────────
# Get the directory where this config file is located (experiments/tcga/)
_CONFIG_DIR = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.dirname(os.path.dirname(_CONFIG_DIR))  # Two levels up to repo root

# GCP credentials: use env var if set, otherwise look in infra/secrets/
_DEFAULT_CREDS = os.path.join(_REPO_ROOT, "infra", "secrets", "gcp-credentials.json")
CREDENTIALS_PATH = os.environ.get("GOOGLE_APPLICATION_CREDENTIALS", _DEFAULT_CREDS)

PROJECT_ID = "kernel-o6"
GCS_BUCKET = "kernel-o6"
GCS_PREFIX = "tcga_cancer_boundary"

# Data and figure directories relative to this script's location
DATA_DIR = _CONFIG_DIR  # Same directory as the scripts
FIGURE_DIR = os.path.join(_CONFIG_DIR, "figures")

# ── Gene Groups ────────────────────────────────────────────────────────────
# Original gene names (for BigQuery queries)
TARGET_GENES_RAW = [
    "CD274", "PDCD1LG2", "PDCD1",
    "HLA-A", "HLA-B", "HLA-C", "B2M",
    "GJA1", "GJB2", "GJA5", "GJB6",
    "ARNTL", "CLOCK", "PER1", "PER2", "CRY1", "CRY2",
    "MYC", "TP53", "CDH1", "VIM",
]

# Column aliases (hyphens replaced for SQL/pandas compatibility)
CHECKPOINT = ["CD274", "PDCD1LG2", "PDCD1"]
MHC_I = ["HLA_A", "HLA_B", "HLA_C", "B2M"]
GAP_JUNCTION = ["GJA1", "GJB2", "GJA5", "GJB6"]
CIRCADIAN = ["ARNTL", "CLOCK", "PER1", "PER2", "CRY1", "CRY2"]
DIFFERENTIATION = ["CDH1", "VIM", "MYC"]

ALL_GENE_COLS = CHECKPOINT + MHC_I + GAP_JUNCTION + CIRCADIAN + DIFFERENTIATION + ["TP53"]

# ── Cancer Types ───────────────────────────────────────────────────────────
ALL_PROJECTS = [
    "TCGA-SKCM",  # Melanoma
    "TCGA-LUAD",  # Lung adenocarcinoma
    "TCGA-BRCA",  # Breast invasive carcinoma
    "TCGA-COAD",  # Colon adenocarcinoma
    "TCGA-HNSC",  # Head and neck squamous cell carcinoma
    "TCGA-LUSC",  # Lung squamous cell carcinoma
]

PROJECT_LABELS = {
    "TCGA-SKCM": "Melanoma",
    "TCGA-LUAD": "Lung Adeno",
    "TCGA-BRCA": "Breast",
    "TCGA-COAD": "Colon",
    "TCGA-HNSC": "Head & Neck",
    "TCGA-LUSC": "Lung Squam.",
}


# ── Utility Functions ──────────────────────────────────────────────────────
def safe_alias(gene_name):
    """Convert gene name to valid SQL/column alias."""
    return gene_name.replace("-", "_")


def log_transform(series):
    """Log2(TPM + 1) transform."""
    return np.log2(series + 1)


def compute_circadian_cv(df_subset):
    """
    Per-sample circadian coherence proxy: coefficient of variation
    across 6 core clock genes. Lower CV = more coherent.
    """
    circ_cols = [c for c in CIRCADIAN if c in df_subset.columns]
    circ_vals = df_subset[circ_cols].apply(log_transform)
    return circ_vals.std(axis=1) / circ_vals.mean(axis=1)


def classify_boundary_failure(df, project_col="project_short_name"):
    """
    Classify tumors into boundary-failure subtypes per cancer type.
    Uses within-project medians for thresholds.

    Returns a Series with labels: Active_Masking, Decoherence, Mixed
    """
    labels = pd.Series("Mixed", index=df.index)

    for proj in df[project_col].unique():
        mask = df[project_col] == proj
        sub = df.loc[mask]

        cd274_med = log_transform(sub["CD274"]).median()
        arntl_med = log_transform(sub["ARNTL"]).median()
        per1_med = log_transform(sub["PER1"]).median()
        b2m_med = log_transform(sub["B2M"]).median()

        cd274_log = log_transform(sub["CD274"])
        arntl_log = log_transform(sub["ARNTL"])
        per1_log = log_transform(sub["PER1"])
        b2m_log = log_transform(sub["B2M"])

        # Active Masking: high PD-L1 + high BMAL1 + low PER1
        am = (cd274_log > cd274_med) & (arntl_log > arntl_med) & (per1_log < per1_med)
        # Decoherence: low PD-L1 + low MHC-I (B2M)
        dc = (cd274_log < cd274_med) & (b2m_log < b2m_med)

        labels.loc[mask & am] = "Active_Masking"
        labels.loc[mask & dc] = "Decoherence"

    return labels


def normalize_stage(stage_str):
    """Normalize TCGA stage labels to I/II/III/IV."""
    if pd.isna(stage_str):
        return None
    s = str(stage_str).lower().strip()
    if "not reported" in s or "unknown" in s:
        return None
    if "iv" in s:
        return "IV"
    if "iii" in s:
        return "III"
    if "ii" in s:
        return "II"
    if "i" in s:
        return "I"
    return None


def get_bq_client():
    """Get authenticated BigQuery client."""
    os.environ["GOOGLE_APPLICATION_CREDENTIALS"] = CREDENTIALS_PATH
    from google.cloud import bigquery
    return bigquery.Client(project=PROJECT_ID)


def get_gcs_bucket():
    """Get GCS bucket handle."""
    os.environ["GOOGLE_APPLICATION_CREDENTIALS"] = CREDENTIALS_PATH
    from google.cloud import storage
    return storage.Client(project=PROJECT_ID).bucket(GCS_BUCKET)


def upload_to_gcs(local_path, gcs_subpath):
    """Upload a local file to GCS under the project prefix."""
    bucket = get_gcs_bucket()
    blob = bucket.blob(f"{GCS_PREFIX}/{gcs_subpath}")
    blob.upload_from_filename(local_path)
    print(f"  -> gs://{GCS_BUCKET}/{GCS_PREFIX}/{gcs_subpath}")


def setup_plotting():
    """Set matplotlib defaults for publication figures."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "figure.dpi": 100,
        "savefig.dpi": 300,
        "font.size": 10,
        "axes.titlesize": 11,
        "axes.labelsize": 10,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.bbox": "tight",
    })
    os.makedirs(FIGURE_DIR, exist_ok=True)
    return plt
