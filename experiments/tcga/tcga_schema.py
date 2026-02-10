"""Check schema and sample data from TCGA RNA-seq table."""
import os
os.environ["GOOGLE_APPLICATION_CREDENTIALS"] = r"C:\Users\Moses\commercials\infra\secrets\gcp-credentials.json"

from google.cloud import bigquery

client = bigquery.Client(project="kernel-o6")

# 1. Get column names
print("=== TABLE SCHEMA ===")
query_schema = """
SELECT column_name, data_type
FROM `isb-cgc-bq.TCGA.INFORMATION_SCHEMA.COLUMNS`
WHERE table_name = 'RNAseq_hg38_gdc_current'
ORDER BY ordinal_position
"""
for row in client.query(query_schema).result():
    print(f"  {row.column_name:30s} {row.data_type}")

# 2. Sample rows
print("\n=== SAMPLE DATA (5 rows, TCGA-SKCM) ===")
query_sample = """
SELECT *
FROM `isb-cgc-bq.TCGA.RNAseq_hg38_gdc_current`
WHERE project_short_name = 'TCGA-SKCM'
  AND gene_name = 'CD274'
LIMIT 5
"""
df = client.query(query_sample).to_dataframe()
print(df.to_string())

# 3. Count samples per project
print("\n=== SAMPLE COUNTS ===")
query_counts = """
SELECT project_short_name, COUNT(DISTINCT case_barcode) AS n_cases
FROM `isb-cgc-bq.TCGA.RNAseq_hg38_gdc_current`
WHERE project_short_name IN ('TCGA-SKCM', 'TCGA-LUAD')
GROUP BY project_short_name
"""
for row in client.query(query_counts).result():
    print(f"  {row.project_short_name}: {row.n_cases} cases")
