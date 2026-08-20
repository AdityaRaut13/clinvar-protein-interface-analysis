#!/usr/bin/env python

import gzip
import json
import os
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
import requests
from tqdm import tqdm


ROOT = Path(os.environ["IF_PROJECT"])

MANIFEST_FILE = (
    ROOT
    / "data/processed/assembly_discovery"
    / "target_complex_assembly_manifest.tsv.gz"
)
OUTPUT_DIR = ROOT / "structures/assemblies"
RESULT_DIR = ROOT / "data/processed/assembly_download"

LIMIT = int(os.environ.get("ASSEMBLY_LIMIT", "0"))
MAX_WORKERS = int(os.environ.get("ASSEMBLY_WORKERS", "8"))
SUFFIX = f"_test{LIMIT}" if LIMIT else ""

STATUS_OUTPUT = RESULT_DIR / f"assembly_download_status{SUFFIX}.tsv.gz"
SUMMARY_OUTPUT = RESULT_DIR / f"assembly_download_summary{SUFFIX}.json"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
RESULT_DIR.mkdir(parents=True, exist_ok=True)

thread_local = threading.local()


def get_session():
    if not hasattr(thread_local, "session"):
        session = requests.Session()
        session.headers.update({
            "User-Agent": "IITK-ClinVar-interface/1.0"
        })
        thread_local.session = session

    return thread_local.session


def validate_gzip(path):
    try:
        if not path.exists() or path.stat().st_size == 0:
            return False

        decompressed_bytes = 0

        with gzip.open(path, "rb") as handle:
            while True:
                chunk = handle.read(1024 * 1024)

                if not chunk:
                    break

                decompressed_bytes += len(chunk)

        return decompressed_bytes > 0

    except (OSError, EOFError):
        return False


def result(record, status, error=""):
    path = OUTPUT_DIR / record.filename

    return {
        "pdb_id": record.pdb_id,
        "assembly_id": record.assembly_id,
        "filename": record.filename,
        "status": status,
        "size_bytes": (
            path.stat().st_size
            if path.exists()
            else 0
        ),
        "error": error,
    }


def download(record):
    destination = OUTPUT_DIR / record.filename
    temporary = OUTPUT_DIR / f"{record.filename}.part"

    if validate_gzip(destination):
        return result(record, "cached")

    destination.unlink(missing_ok=True)
    temporary.unlink(missing_ok=True)

    last_error = ""

    for attempt in range(6):
        try:
            with get_session().get(
                record.download_url,
                stream=True,
                timeout=300,
            ) as response:
                if response.status_code == 404:
                    return result(
                        record,
                        "not_found",
                        "HTTP 404",
                    )

                response.raise_for_status()

                with open(temporary, "wb") as handle:
                    for chunk in response.iter_content(
                        chunk_size=1024 * 1024
                    ):
                        if chunk:
                            handle.write(chunk)

            if not validate_gzip(temporary):
                raise OSError(
                    "Downloaded file failed gzip validation"
                )

            temporary.replace(destination)
            return result(record, "downloaded")

        except (requests.RequestException, OSError) as error:
            last_error = str(error)
            temporary.unlink(missing_ok=True)
            time.sleep(2 ** attempt)

    return result(record, "error", last_error)


manifest = pd.read_csv(
    MANIFEST_FILE,
    sep="\t",
    dtype={
        "pdb_id": str,
        "assembly_id": str,
    },
)

manifest = (
    manifest.drop_duplicates(
        ["pdb_id", "assembly_id"]
    )
    .sort_values(["pdb_id", "assembly_id"])
)

if LIMIT:
    manifest = manifest.head(LIMIT)

print(f"Assemblies to process: {len(manifest):,}")
print(f"Download workers: {MAX_WORKERS}")

records = list(manifest.itertuples(index=False))
rows = []

with ThreadPoolExecutor(
    max_workers=MAX_WORKERS
) as executor:
    futures = {
        executor.submit(download, record): record
        for record in records
    }

    for future in tqdm(
        as_completed(futures),
        total=len(futures),
        desc="Downloading biological assemblies",
    ):
        record = futures[future]

        try:
            rows.append(future.result())
        except Exception as error:
            rows.append(
                result(record, "worker_error", str(error))
            )

status = pd.DataFrame(rows).sort_values(
    ["status", "pdb_id", "assembly_id"]
)

status.to_csv(
    STATUS_OUTPUT,
    sep="\t",
    index=False,
    compression="gzip",
)

successful = status["status"].isin({
    "downloaded",
    "cached",
})

summary = {
    "requested_assemblies": int(len(manifest)),
    "successful_assemblies": int(successful.sum()),
    "failed_assemblies": int((~successful).sum()),
    "downloaded": int(
        status["status"].eq("downloaded").sum()
    ),
    "cached": int(
        status["status"].eq("cached").sum()
    ),
    "not_found": int(
        status["status"].eq("not_found").sum()
    ),
    "errors": int(
        status["status"].isin({
            "error",
            "worker_error",
        }).sum()
    ),
    "total_size_bytes": int(
        status.loc[successful, "size_bytes"].sum()
    ),
    "total_size_gib": round(
        status.loc[successful, "size_bytes"].sum()
        / (1024 ** 3),
        3,
    ),
}

with open(SUMMARY_OUTPUT, "w") as handle:
    json.dump(summary, handle, indent=2)

print()
print(json.dumps(summary, indent=2))
print()
print("Status:", STATUS_OUTPUT)