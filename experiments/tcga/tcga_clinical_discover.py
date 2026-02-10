"""Discover clinical table schema and sample data."""
import os, sys
sys.path.insert(0, os.path.dirname(__file__))
from tcga_config import get_bq_client

client = get_bq_client()

print("=== CLINICAL TABLE COLUMNS ===")
q = """
SELECT column_name, data_type
FROM `isb-cgc-bq.TCGA.INFORMATION_SCHEMA.COLUMNS`
WHERE table_name = 'clinical_gdc_current'
ORDER BY ordinal_position
"""
for row in client.query(q).result():
    print(f"  {row.column_name:40s} {row.data_type}")

print("\n=== SAMPLE ROW ===")
q2 = """
SELECT * FROM `isb-cgc-bq.TCGA.clinical_gdc_current`
WHERE proj__project_id = 'TCGA-SKCM'
LIMIT 2
"""
df = client.query(q2).to_dataframe()
for col in df.columns:
    print(f"  {col:40s} = {df[col].iloc[0]}")
