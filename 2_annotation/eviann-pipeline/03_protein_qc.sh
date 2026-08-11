#!/usr/bin/env bash

# protein-seq for eviann
# EviAnn's -p input = single FASTA 

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_config.sh"

PROTEINS_OUT="${OUT_PROTQC}/proteins.clean.faa"
MANIFEST="${OUT_PROTQC}/proteins.clean.manifest.tsv"
REPORT="${OUT_PROTQC}/proteins.clean.report.txt"

mkdir -p "${OUT_PROTQC}"

python3 - "${PROT_DIR}" "${PROTEINS_OUT}" "${MANIFEST}" "${REPORT}" <<'PY'
from __future__ import annotations

import gzip
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path


prot_dir = Path(sys.argv[1]).resolve()
proteins_out = Path(sys.argv[2]).resolve()
manifest = Path(sys.argv[3]).resolve()
report = Path(sys.argv[4]).resolve()

suffixes = (".faa", ".fa", ".fasta", ".fsa", ".faa.gz", ".fa.gz", ".fasta.gz", ".fsa.gz")
valid_gap_or_stop = set("-.*")
valid_protein = set("ABCDEFGHIKLMNPQRSTUVWYZXJOUB")


def fasta_inputs(directory: Path) -> list[Path]:
    files = [
        path
        for path in directory.iterdir()
        if path.is_file() and not path.name.startswith(".") and path.name.lower().endswith(suffixes)
    ]
    return sorted(files, key=lambda path: path.name)


def open_text(path: Path):
    if path.name.lower().endswith(".gz"):
        return gzip.open(path, "rt")
    return path.open()


def parse_fasta(path: Path):
    header = None
    seq_chunks: list[str] = []
    with open_text(path) as handle:
        for line_no, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks)
                header = line[1:].strip()
                seq_chunks = []
            elif header is None:
                raise ValueError(f"{path}: sequence data before first FASTA header at line {line_no}")
            else:
                seq_chunks.append(line)
    if header is not None:
        yield header, "".join(seq_chunks)


def clean_id(raw_id: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.:|+-]", "_", raw_id)
    cleaned = cleaned.strip("_")
    return cleaned or "protein"


def clean_sequence(raw_sequence: str) -> tuple[str, Counter[str]]:
    replacements: Counter[str] = Counter()
    cleaned: list[str] = []
    for char in raw_sequence.upper():
        if char.isspace() or char.isdigit():
            replacements[char] += 1
            continue
        if char in valid_gap_or_stop:
            replacements[char] += 1
            continue
        if char in valid_protein:
            cleaned.append(char)
        elif char.isalpha():
            cleaned.append("X")
            replacements[char] += 1
        else:
            replacements[char] += 1
    return "".join(cleaned), replacements


inputs = fasta_inputs(prot_dir)
if not prot_dir.is_dir():
    raise SystemExit(f"ERROR: protein directory does not exist: {prot_dir}")
if not inputs:
    raise SystemExit(f"ERROR: no FASTA inputs found in {prot_dir}")

seen_ids: Counter[str] = Counter()
summary_rows: list[dict[str, object]] = []
global_replacements: Counter[str] = Counter()
total_in = 0
total_written = 0
total_empty = 0
total_renamed = 0
total_bases = 0
length_min: int | None = None
length_max = 0
source_counts: defaultdict[str, int] = defaultdict(int)

with proteins_out.open("w") as out:
    for input_path in inputs:
        records_in = 0
        records_written = 0
        skipped_empty = 0
        renamed = 0
        replacements: Counter[str] = Counter()

        for header, raw_sequence in parse_fasta(input_path):
            records_in += 1
            total_in += 1

            first_token, *rest = header.split(maxsplit=1)
            description = rest[0] if rest else ""
            base_id = clean_id(first_token)
            seen_ids[base_id] += 1
            output_id = base_id
            if seen_ids[base_id] > 1:
                source_counts[input_path.stem] += 1
                output_id = f"{base_id}|{input_path.stem}|dup{seen_ids[base_id]}"
                renamed += 1
                total_renamed += 1

            sequence, record_replacements = clean_sequence(raw_sequence)
            replacements.update(record_replacements)
            global_replacements.update(record_replacements)
            if not sequence:
                skipped_empty += 1
                total_empty += 1
                continue

            seq_len = len(sequence)
            total_bases += seq_len
            length_min = seq_len if length_min is None else min(length_min, seq_len)
            length_max = max(length_max, seq_len)
            total_written += 1
            records_written += 1

            clean_description = " ".join(description.replace("\t", " ").split())
            out.write(f">{output_id}")
            if clean_description:
                out.write(f" {clean_description}")
            out.write("\n")
            for start in range(0, seq_len, 60):
                out.write(sequence[start:start + 60] + "\n")

        summary_rows.append(
            {
                "file": str(input_path),
                "records_in": records_in,
                "records_written": records_written,
                "skipped_empty": skipped_empty,
                "duplicate_ids_renamed": renamed,
                "chars_removed_or_replaced": sum(replacements.values()),
            }
        )

with manifest.open("w") as handle:
    handle.write(
        "file\trecords_in\trecords_written\tskipped_empty\t"
        "duplicate_ids_renamed\tchars_removed_or_replaced\n"
    )
    for row in summary_rows:
        handle.write(
            "{file}\t{records_in}\t{records_written}\t{skipped_empty}\t"
            "{duplicate_ids_renamed}\t{chars_removed_or_replaced}\n".format(**row)
        )

with report.open("w") as handle:
    handle.write("EviAnn protein evidence preparation report\n")
    handle.write(f"Input directory: {prot_dir}\n")
    handle.write(f"Output FASTA: {proteins_out}\n")
    handle.write(f"Manifest: {manifest}\n")
    handle.write(f"Input FASTA files: {len(inputs)}\n")
    handle.write(f"Records read: {total_in}\n")
    handle.write(f"Records written: {total_written}\n")
    handle.write(f"Empty records skipped: {total_empty}\n")
    handle.write(f"Duplicate FASTA IDs renamed: {total_renamed}\n")
    handle.write(f"Total amino-acid residues: {total_bases}\n")
    handle.write(f"Minimum protein length: {length_min or 0}\n")
    handle.write(f"Maximum protein length: {length_max}\n")
    if global_replacements:
        handle.write("Characters removed/replaced:\n")
        for char, count in sorted(global_replacements.items(), key=lambda item: repr(item[0])):
            label = {" ": "space", "\t": "tab", "\n": "newline"}.get(char, char)
            handle.write(f"  {label}: {count}\n")

print(f"Wrote {total_written:,} protein records to {proteins_out}")
print(f"Manifest: {manifest}")
print(f"Report: {report}")
PY

if [[ ! -s "${PROTEINS_OUT}" ]]; then
    echo "ERROR: protein output was not created: ${PROTEINS_OUT}" >&2
    exit 1
fi

echo "[$(date)] Protein evidence ready for EviAnn -p: ${PROTEINS_OUT}"
