"""Discover available TCGA RNA-seq tables in ISB-CGC BigQuery."""
import os
os.environ["GOOGLE_APPLICATION_CREDENTIALS"] = r"C:\Users\Moses\commercials\infra\secrets\gcp-credentials.json"

from google.cloud import bigquery

client = bigquery.Client(project="kernel-o6")

# Check what TCGA RNA-seq tables exist
query = """
SELECT table_name
FROM `isb-cgc-bq.TCGA.INFORMATION_SCHEMA.TABLES`
WHERE table_name LIKE '%RNAseq%'
ORDER BY table_name
"""

print("Querying ISB-CGC BigQuery for TCGA RNA-seq tables...")
results = client.query(query).result()
for row in results:
    print(f"  {row.table_name}")
