"""Discover available TCGA RNA-seq tables in ISB-CGC BigQuery."""
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
from tcga_config import get_bq_client  # noqa: E402

client = get_bq_client()

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
