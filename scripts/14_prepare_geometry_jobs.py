#!/usr/bin/env python3

import json
import os
from pathlib import Path

import duckdb

PROJECT = Path(os.environ["IF_PROJECT"])

TASKS = (
    PROJECT / "data/processed/variant_pair_tasks/"
    "variant_chain_pair_tasks.parquet"
)

OUTDIR = PROJECT / "data/processed/variant_geometry"
JOBS = OUTDIR / "geometry_jobs"
SUMMARY = OUTDIR / "geometry_job_summary.json"

TEMP = Path(os.environ["IF_TEMP"])

OUTDIR.mkdir(parents=True, exist_ok=True)
TEMP.mkdir(parents=True, exist_ok=True)

if JOBS.exists() and any(JOBS.rglob("*.parquet")):
    raise FileExistsError(f"Output already exists: {JOBS}")

con = duckdb.connect()
con.execute("SET threads = 8")
con.execute("SET memory_limit = '24GB'")
con.execute(f"SET temp_directory = '{TEMP}'")
con.execute("SET preserve_insertion_order = false")
con.execute("PRAGMA enable_progress_bar")

print("Grouping variant positions by chain pair", flush=True)

con.execute(
    f"""
    COPY (
        WITH grouped AS (
            SELECT
                pdb_id,
                assembly_id,
                target_label_asym_id,
                partner_label_asym_id,
                list(
                    DISTINCT pdb_label_seq_id
                    ORDER BY pdb_label_seq_id
                ) AS pdb_label_seq_ids
            FROM read_parquet('{TASKS}')
            GROUP BY
                pdb_id,
                assembly_id,
                target_label_asym_id,
                partner_label_asym_id
        )

        SELECT
            hash(pdb_id, assembly_id) % 128 AS task_bucket,
            *
        FROM grouped
    )
    TO '{JOBS}'
    (
        FORMAT PARQUET,
        COMPRESSION ZSTD,
        PARTITION_BY (task_bucket)
    )
    """
)

job_glob = JOBS / "**/*.parquet"

row = con.execute(
    f"""
    SELECT
        count(*) AS chain_pair_jobs,
        count(DISTINCT pdb_id || ':' || assembly_id) AS assemblies,
        sum(list_count(pdb_label_seq_ids)) AS residue_pair_tasks,
        max(list_count(pdb_label_seq_ids)) AS maximum_positions_per_pair
    FROM read_parquet('{job_glob}')
    """
).fetchone()

summary = {
    "chain_pair_jobs": int(row[0]),
    "assemblies": int(row[1]),
    "residue_pair_tasks": int(row[2]),
    "maximum_positions_per_chain_pair": int(row[3]),
    "task_buckets": 128,
}

with open(SUMMARY, "w") as handle:
    json.dump(summary, handle, indent=2)

print(json.dumps(summary, indent=2))
print("\nGeometry jobs:", JOBS)
print("Summary:", SUMMARY)

con.close()