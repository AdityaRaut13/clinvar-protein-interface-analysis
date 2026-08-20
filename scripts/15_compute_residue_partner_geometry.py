#!/usr/bin/env python3

import gzip
import json
import os
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import duckdb
import gemmi
import numpy as np
import pyarrow as pa
import pyarrow.dataset as ds
import pyarrow.parquet as pq
from scipy.spatial import cKDTree
from tqdm import tqdm


PROJECT = Path(os.environ["IF_PROJECT"])

JOBS = PROJECT / "data/processed/variant_geometry/geometry_jobs"
ASSEMBLIES = PROJECT / "structures/assemblies"

OUTDIR = PROJECT / "data/processed/variant_geometry"
GEOMETRY = OUTDIR / "residue_partner_geometry_union10A5A"
STATUS = OUTDIR / "assembly_geometry_status_union10A5A"
SUMMARY = OUTDIR / "residue_partner_geometry_summary_union10A5A.json"

WORKERS = int(os.environ.get("IF_GEOMETRY_WORKERS", "8"))
MAX_BUCKETS = int(os.environ.get("IF_GEOMETRY_MAX_BUCKETS", "0"))

BACKBONE_ATOMS = {"N", "CA", "C", "O"}

GEOMETRY_SCHEMA = pa.schema([
    ("pdb_id", pa.string()),
    ("assembly_id", pa.string()),
    ("target_label_asym_id", pa.string()),
    ("partner_label_asym_id", pa.string()),
    ("pdb_label_seq_id", pa.int32()),
    ("pdb_residue_name", pa.string()),
    ("n_target_heavy_atoms", pa.int16()),
    ("n_target_backbone_atoms", pa.int16()),
    ("min_backbone_distance_to_partner", pa.float32()),
    ("min_heavy_atom_distance_to_partner", pa.float32()),
    ("residue_pinder_interface_10A", pa.bool_()),
    ("residue_heavy_interface_5A", pa.bool_()),
    ("geometry_status", pa.string()),
])

STATUS_SCHEMA = pa.schema([
    ("pdb_id", pa.string()),
    ("assembly_id", pa.string()),
    ("status", pa.string()),
    ("n_chain_pair_jobs", pa.int64()),
    ("n_residue_tasks", pa.int64()),
    ("n_geometry_rows", pa.int64()),
    ("error", pa.string()),
])


def values(block, tag):
    return list(block.find_values(tag))


def optional(items, index, default=""):
    if index >= len(items):
        return default

    value = str(items[index]).strip()

    if value in {"", ".", "?"}:
        return default

    return value


def clean_id(value):
    value = str(value).strip()
    return "" if value in {"", ".", "?"} else value


def model_sort_key(value):
    try:
        return 0, int(value)
    except ValueError:
        return 1, value


def parse_needed_chains(path, needed_chains):
    with gzip.open(
        path,
        "rt",
        encoding="utf-8",
        errors="replace",
    ) as handle:
        document = gemmi.cif.read_string(handle.read())

    block = document.sole_block()

    label_asym = values(block, "_atom_site.label_asym_id")
    label_seq = values(block, "_atom_site.label_seq_id")
    residue_names = values(block, "_atom_site.label_comp_id")
    elements = values(block, "_atom_site.type_symbol")
    atom_names = values(block, "_atom_site.label_atom_id")
    x_values = values(block, "_atom_site.Cartn_x")
    y_values = values(block, "_atom_site.Cartn_y")
    z_values = values(block, "_atom_site.Cartn_z")
    models = values(block, "_atom_site.pdbx_PDB_model_num")

    if not label_asym:
        raise ValueError("Missing atom_site records")

    if not models:
        models = ["1"] * len(label_asym)

    model_used = sorted(set(models), key=model_sort_key)[0]

    chains = defaultdict(lambda: {
        "heavy": [],
        "backbone": [],
        "residues": defaultdict(lambda: {
            "name": "",
            "heavy": [],
            "backbone": [],
        }),
    })

    for index, raw_label in enumerate(label_asym):
        if optional(models, index, "1") != model_used:
            continue

        label_id = clean_id(raw_label)

        if label_id not in needed_chains:
            continue

        element = optional(elements, index).upper()

        if element in {"H", "D"}:
            continue

        try:
            coordinate = (
                float(x_values[index]),
                float(y_values[index]),
                float(z_values[index]),
            )
        except (ValueError, IndexError):
            continue

        atom_name = optional(atom_names, index).upper()
        chain = chains[label_id]

        chain["heavy"].append(coordinate)

        if atom_name in BACKBONE_ATOMS:
            chain["backbone"].append(coordinate)

        seq_value = optional(label_seq, index)

        try:
            seq_id = int(seq_value)
        except ValueError:
            continue

        residue = chain["residues"][seq_id]

        if not residue["name"]:
            residue["name"] = optional(residue_names, index)

        residue["heavy"].append(coordinate)

        if atom_name in BACKBONE_ATOMS:
            residue["backbone"].append(coordinate)

    for chain in chains.values():
        chain["heavy"] = np.asarray(
            chain["heavy"], dtype=np.float32
        ).reshape(-1, 3)

        chain["backbone"] = np.asarray(
            chain["backbone"], dtype=np.float32
        ).reshape(-1, 3)

        for residue in chain["residues"].values():
            residue["heavy"] = np.asarray(
                residue["heavy"], dtype=np.float32
            ).reshape(-1, 3)

            residue["backbone"] = np.asarray(
                residue["backbone"], dtype=np.float32
            ).reshape(-1, 3)

    return dict(chains)


def minimum_distance(coordinates, tree):
    if coordinates.shape[0] == 0 or tree is None:
        return None

    distances, _ = tree.query(coordinates, k=1)
    return float(np.min(distances))


def process_bucket(bucket_number):
    input_dir = JOBS / f"task_bucket={bucket_number}"
    geometry_output = GEOMETRY / f"part-{bucket_number:03d}.parquet"
    status_output = STATUS / f"part-{bucket_number:03d}.parquet"

    table = ds.dataset(input_dir, format="parquet").to_table(
        columns=[
            "pdb_id",
            "assembly_id",
            "target_label_asym_id",
            "partner_label_asym_id",
            "pdb_label_seq_ids",
        ]
    )
    frame = table.to_pandas()

    geometry_records = []
    status_records = []

    for (pdb_id, assembly_id), group in frame.groupby(
        ["pdb_id", "assembly_id"],
        sort=False,
    ):
        pdb_id = str(pdb_id).lower()
        assembly_id = str(assembly_id)

        path = ASSEMBLIES / (
            f"{pdb_id}-assembly{assembly_id}.cif.gz"
        )

        n_tasks = int(
            sum(len(x) for x in group["pdb_label_seq_ids"])
        )

        needed_chains = set(group["target_label_asym_id"])
        needed_chains.update(group["partner_label_asym_id"])

        try:
            if not path.exists():
                raise FileNotFoundError(path)

            chains = parse_needed_chains(path, needed_chains)
            tree_cache = {}

            def get_trees(chain_id):
                if chain_id in tree_cache:
                    return tree_cache[chain_id]

                chain = chains.get(chain_id)

                if chain is None:
                    result = (None, None)
                else:
                    heavy_tree = (
                        cKDTree(chain["heavy"])
                        if len(chain["heavy"]) else None
                    )
                    backbone_tree = (
                        cKDTree(chain["backbone"])
                        if len(chain["backbone"]) else None
                    )
                    result = heavy_tree, backbone_tree

                tree_cache[chain_id] = result
                return result

            assembly_rows = 0
            missing_rows = 0

            for row in group.itertuples(index=False):
                target = chains.get(row.target_label_asym_id)
                partner = chains.get(row.partner_label_asym_id)

                partner_heavy_tree, partner_backbone_tree = (
                    get_trees(row.partner_label_asym_id)
                )

                for raw_position in row.pdb_label_seq_ids:
                    position = int(raw_position)
                    residue = (
                        target["residues"].get(position)
                        if target is not None else None
                    )

                    heavy_distance = None
                    backbone_distance = None
                    residue_name = None
                    n_heavy = 0
                    n_backbone = 0

                    if target is None:
                        geometry_status = "target_chain_missing"
                    elif partner is None:
                        geometry_status = "partner_chain_missing"
                    elif residue is None:
                        geometry_status = "target_residue_missing"
                    else:
                        residue_name = residue["name"]
                        n_heavy = len(residue["heavy"])
                        n_backbone = len(residue["backbone"])

                        heavy_distance = minimum_distance(
                            residue["heavy"],
                            partner_heavy_tree,
                        )
                        backbone_distance = minimum_distance(
                            residue["backbone"],
                            partner_backbone_tree,
                        )

                        if heavy_distance is None:
                            geometry_status = "heavy_geometry_missing"
                        elif backbone_distance is None:
                            geometry_status = "backbone_geometry_missing"
                        else:
                            geometry_status = "ok"

                    if geometry_status != "ok":
                        missing_rows += 1

                    geometry_records.append({
                        "pdb_id": pdb_id,
                        "assembly_id": assembly_id,
                        "target_label_asym_id":
                            row.target_label_asym_id,
                        "partner_label_asym_id":
                            row.partner_label_asym_id,
                        "pdb_label_seq_id": position,
                        "pdb_residue_name": residue_name,
                        "n_target_heavy_atoms": n_heavy,
                        "n_target_backbone_atoms": n_backbone,
                        "min_backbone_distance_to_partner":
                            backbone_distance,
                        "min_heavy_atom_distance_to_partner":
                            heavy_distance,
                        "residue_pinder_interface_10A":
                            (
                                backbone_distance <= 10.0
                                if backbone_distance is not None
                                else None
                            ),
                        "residue_heavy_interface_5A":
                            (
                                heavy_distance <= 5.0
                                if heavy_distance is not None
                                else None
                            ),
                        "geometry_status": geometry_status,
                    })

                    assembly_rows += 1

            status = (
                "ok"
                if missing_rows == 0
                else "ok_with_missing_geometry"
            )

            status_records.append({
                "pdb_id": pdb_id,
                "assembly_id": assembly_id,
                "status": status,
                "n_chain_pair_jobs": len(group),
                "n_residue_tasks": n_tasks,
                "n_geometry_rows": assembly_rows,
                "error": "",
            })

        except Exception as error:
            status_records.append({
                "pdb_id": pdb_id,
                "assembly_id": assembly_id,
                "status": "error",
                "n_chain_pair_jobs": len(group),
                "n_residue_tasks": n_tasks,
                "n_geometry_rows": 0,
                "error": str(error)[:1000],
            })

    geometry_table = pa.Table.from_pylist(
        geometry_records,
        schema=GEOMETRY_SCHEMA,
    )
    status_table = pa.Table.from_pylist(
        status_records,
        schema=STATUS_SCHEMA,
    )

    geometry_tmp = geometry_output.with_suffix(".tmp.parquet")
    status_tmp = status_output.with_suffix(".tmp.parquet")

    pq.write_table(
        geometry_table,
        geometry_tmp,
        compression="zstd",
        row_group_size=100000,
    )
    pq.write_table(
        status_table,
        status_tmp,
        compression="zstd",
    )

    geometry_tmp.replace(geometry_output)
    status_tmp.replace(status_output)

    return bucket_number


def build_summary():
    geometry_glob = GEOMETRY / "*.parquet"
    status_glob = STATUS / "*.parquet"
    jobs_glob = JOBS / "**/*.parquet"

    con = duckdb.connect()
    con.execute("SET threads = 8")

    expected = con.execute(
        f"""
        SELECT sum(list_count(pdb_label_seq_ids))
        FROM read_parquet('{jobs_glob}')
        """
    ).fetchone()[0]

    completed = len(list(GEOMETRY.glob("part-*.parquet")))

    if completed:
        row = con.execute(
            f"""
            SELECT
                count(*),
                count(*) FILTER (
                    WHERE geometry_status = 'ok'
                ),
                count(*) FILTER (
                    WHERE residue_pinder_interface_10A
                ),
                count(*) FILTER (
                    WHERE residue_heavy_interface_5A
                ),
                count(*) FILTER (
                    WHERE residue_pinder_interface_10A
                      AND residue_heavy_interface_5A
                ),
                count(*) FILTER (
                    WHERE residue_pinder_interface_10A
                      AND NOT residue_heavy_interface_5A
                ),
                count(*) FILTER (
                    WHERE NOT residue_pinder_interface_10A
                      AND residue_heavy_interface_5A
                )
            FROM read_parquet('{geometry_glob}')
            """
        ).fetchone()

        geometry_status = dict(
            con.execute(
                f"""
                SELECT geometry_status, count(*)
                FROM read_parquet('{geometry_glob}')
                GROUP BY geometry_status
                """
            ).fetchall()
        )
    else:
        row = (0,) * 7
        geometry_status = {}

    assembly_status = {}

    if list(STATUS.glob("part-*.parquet")):
        assembly_status = dict(
            con.execute(
                f"""
                SELECT status, count(*)
                FROM read_parquet('{status_glob}')
                GROUP BY status
                """
            ).fetchall()
        )

    summary = {
        "expected_residue_pair_tasks": int(expected),
        "completed_task_buckets": completed,
        "total_task_buckets": 128,
        "geometry_rows": int(row[0]),
        "geometry_ok_rows": int(row[1]),
        "pinder_interface_10A_rows": int(row[2]),
        "heavy_interface_5A_rows": int(row[3]),
        "both_interface_definitions": int(row[4]),
        "pinder_only": int(row[5]),
        "heavy_only": int(row[6]),
        "geometry_status_counts": {
            str(key): int(value)
            for key, value in geometry_status.items()
        },
        "assembly_status_counts": {
            str(key): int(value)
            for key, value in assembly_status.items()
        },
    }

    with open(SUMMARY, "w") as handle:
        json.dump(summary, handle, indent=2)

    print(json.dumps(summary, indent=2))
    con.close()


def main():
    GEOMETRY.mkdir(parents=True, exist_ok=True)
    STATUS.mkdir(parents=True, exist_ok=True)

    bucket_numbers = sorted(
        int(path.name.split("=", 1)[1])
        for path in JOBS.glob("task_bucket=*")
    )

    pending = [
        bucket
        for bucket in bucket_numbers
        if not (
            GEOMETRY / f"part-{bucket:03d}.parquet"
        ).exists()
    ]

    if MAX_BUCKETS > 0:
        pending = pending[:MAX_BUCKETS]

    print("Workers:", WORKERS)
    print("Pending buckets selected:", len(pending))

    if pending:
        with ProcessPoolExecutor(
            max_workers=WORKERS
        ) as executor:
            futures = {
                executor.submit(process_bucket, bucket): bucket
                for bucket in pending
            }

            for future in tqdm(
                as_completed(futures),
                total=len(futures),
                desc="Geometry buckets",
            ):
                future.result()

    build_summary()


if __name__ == "__main__":
    main()