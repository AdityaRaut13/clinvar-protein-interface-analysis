import gzip
import pandas as pd
from Bio import SeqIO
from pathlib import Path
import os

root = Path(os.environ["IF_PROJECT"])
csv_path = root / "data/input/clinvar_missense_protein_ready.csv.gz"
fasta_path = root / "data/input/protein_reference_sequences.fasta.gz"

df = pd.read_csv(
    csv_path,
    usecols=["protein_sequence_id", "protein_id"],
    dtype=str,
)

with gzip.open(fasta_path, "rt") as handle:
    fasta_ids = {record.id for record in SeqIO.parse(handle, "fasta")}

print("FASTA records:", f"{len(fasta_ids):,}")
print("Example FASTA IDs:", sorted(fasta_ids)[:5])

for column in ["protein_sequence_id", "protein_id"]:
    values = set(df[column].dropna())
    matched = values & fasta_ids
    print(
        column,
        "unique =", f"{len(values):,}",
        "matched =", f"{len(matched):,}",
        "coverage =", f"{len(matched) / len(values):.2%}" if values else "NA",
    )