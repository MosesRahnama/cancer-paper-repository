"""Discover diagnosis and follow-up tables for stage + survival time."""
import os, sys
sys.path.insert(0, os.path.dirname(__file__))
from tcga_config import get_bq_client

client = get_bq_client()

# Check what other clinical-related tables exist
print("=== ALL TCGA TABLES ===")
q = """
SELECT table_name
FROM `isb-cgc-bq.TCGA.INFORMATION_SCHEMA.TABLES`
WHERE table_name LIKE '%clinical%' OR table_name LIKE '%diag%' OR table_name LIKE '%follow%'
ORDER BY table_name
"""
for row in client.query(q).result():
    print(f"  {row.table_name}")

# Check for diagnosis table (has staging)
print("\n=== DIAGNOSIS TABLE COLUMNS (if exists) ===")
try:
    q2 = """
    SELECT column_name, data_type
    FROM `isb-cgc-bq.TCGA.INFORMATION_SCHEMA.COLUMNS`
    WHERE table_name = 'clinical_diagnoses_gdc_current'
    ORDER BY ordinal_position
    """
    for row in client.query(q2).result():
        print(f"  {row.column_name:40s} {row.data_type}")
except Exception as e:
    print(f"  Not found: {e}")

# Check for follow-up table
print("\n=== FOLLOW-UP TABLE COLUMNS (if exists) ===")
try:
    q3 = """
    SELECT column_name, data_type
    FROM `isb-cgc-bq.TCGA.INFORMATION_SCHEMA.COLUMNS`
    WHERE table_name = 'clinical_follow_ups_gdc_current'
    ORDER BY ordinal_position
    """
    for row in client.query(q3).result():
        print(f"  {row.column_name:40s} {row.data_type}")
except Exception as e:
    print(f"  Not found: {e}")
