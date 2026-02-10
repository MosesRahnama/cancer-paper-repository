"""Get schemas for diagnosis and follow-up tables."""
import os, sys
sys.path.insert(0, os.path.dirname(__file__))
from tcga_config import get_bq_client

client = get_bq_client()

for table in ["clinical_diagnosis_gdc_current", "clinical_follow_up_gdc_current"]:
    print(f"\n=== {table} ===")
    q = f"""
    SELECT column_name, data_type
    FROM `isb-cgc-bq.TCGA.INFORMATION_SCHEMA.COLUMNS`
    WHERE table_name = '{table}'
    ORDER BY ordinal_position
    """
    for row in client.query(q).result():
        print(f"  {row.column_name:45s} {row.data_type}")

    # Sample row
    print(f"\n  --- Sample (TCGA-SKCM) ---")
    q2 = f"""
    SELECT * FROM `isb-cgc-bq.TCGA.{table}`
    WHERE proj__project_id = 'TCGA-SKCM'
    LIMIT 1
    """
    df = client.query(q2).to_dataframe()
    for col in df.columns:
        val = df[col].iloc[0] if len(df) > 0 else "EMPTY"
        print(f"  {col:45s} = {val}")
