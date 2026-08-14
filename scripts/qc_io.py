"""Input validation and atomic output helpers for metadata QC."""

from __future__ import annotations

import csv
import hashlib
import json
import os
import tempfile
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any


AUDIT_FIELDS = [
    "canonical_id",
    "gcf_accession",
    "gca_accession",
    "pair_status",
    "missing_partner_accession",
    "original_representative_accession",
    "original_representative_status",
    "selected_accession",
    "selected_source",
    "selected_status",
    "selection_reason",
    "organism_name",
    "tax_id",
    "strain",
    "assembly_level",
    "genome_size_bp",
    "gc_percent",
    "number_of_contigs",
    "contig_n50",
    "genome_size_modified_z",
    "gc_percent_modified_z",
    "metadata_decision",
    "pangenome_contiguity_class",
    "synteny_class",
    "sequence_qc_status",
    "taxonomy_qc_status",
    "flags",
]


class MetadataQCError(ValueError):
    """Raised when QC inputs or policy violate the metadata contract."""


def sha256_file(path: Path) -> str:
    """Return the SHA-256 digest of a file."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_policy(path: Path) -> dict[str, Any]:
    """Load and validate the metadata-QC policy."""
    try:
        policy = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise MetadataQCError(f"Cannot read QC policy {path}: {exc}.") from exc
    if not isinstance(policy, dict):
        raise MetadataQCError("QC policy must be a JSON object.")
    required = {
        "policy_id",
        "current_status_values",
        "representative_selection",
        "metadata_contiguity",
        "robust_outlier_detection",
    }
    missing = sorted(required - set(policy))
    if missing:
        raise MetadataQCError(f"QC policy is missing keys: {missing!r}.")
    statuses = policy["current_status_values"]
    if not isinstance(statuses, list) or not statuses:
        raise MetadataQCError("current_status_values must be a non-empty list.")
    contiguity = policy["metadata_contiguity"]
    numeric_keys = (
        "pangenome_max_contigs",
        "pangenome_min_contig_n50_bp",
        "synteny_max_contigs",
        "synteny_min_contig_n50_bp",
    )
    try:
        if any(float(contiguity[key]) <= 0 for key in numeric_keys):
            raise MetadataQCError("All contiguity thresholds must be positive.")
        threshold = float(policy["robust_outlier_detection"]["threshold"])
    except (KeyError, TypeError, ValueError) as exc:
        raise MetadataQCError("QC policy contains invalid numeric thresholds.") from exc
    if threshold <= 0:
        raise MetadataQCError("The robust outlier threshold must be positive.")
    return policy


def load_summary(path: Path) -> list[dict[str, str]]:
    """Load the canonical assembly TSV and reject duplicate identifiers."""
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        fields = set(reader.fieldnames or [])
    required = {
        "canonical_id",
        "representative_accession",
        "gcf_accession",
        "gca_accession",
        "pair_status",
        "missing_partner_accession",
    }
    missing = sorted(required - fields)
    if missing:
        raise MetadataQCError(f"Assembly summary is missing columns: {missing!r}.")
    identifiers = [row["canonical_id"] for row in rows]
    if len(identifiers) != len(set(identifiers)):
        raise MetadataQCError("Assembly summary contains duplicate canonical IDs.")
    if not rows:
        raise MetadataQCError("Assembly summary is empty.")
    return rows


def load_raw_records(path: Path) -> dict[str, dict[str, Any]]:
    """Load raw NCBI JSON Lines records with duplicate-accession protection."""
    records: dict[str, dict[str, Any]] = {}
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                raise MetadataQCError(f"Blank raw record at line {line_number}.")
            try:
                record = json.loads(line)
            except json.JSONDecodeError as exc:
                raise MetadataQCError(
                    f"Invalid JSON at line {line_number}: {exc.msg}."
                ) from exc
            if not isinstance(record, dict):
                raise MetadataQCError(f"Raw record at line {line_number} is not an object.")
            accession = record.get("accession")
            if not isinstance(accession, str) or not accession:
                raise MetadataQCError(f"Missing accession at raw line {line_number}.")
            if accession in records:
                raise MetadataQCError(f"Duplicate raw accession {accession!r}.")
            records[accession] = record
    if not records:
        raise MetadataQCError("Raw metadata input is empty.")
    return records


def write_tsv_atomic(
    rows: Sequence[Mapping[str, Any]], path: Path, *, overwrite: bool
) -> None:
    """Write an audit TSV atomically and refuse silent overwrite."""
    if path.exists() and not overwrite:
        raise FileExistsError(f"Output already exists: {path}. Use --overwrite explicitly.")
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            "w",
            encoding="utf-8",
            newline="",
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            temporary_path = Path(handle.name)
            writer = csv.DictWriter(
                handle,
                fieldnames=AUDIT_FIELDS,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temporary_path, path)
        temporary_path = None
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)


def write_json_atomic(payload: Mapping[str, Any], path: Path, *, overwrite: bool) -> None:
    """Write a JSON report atomically and refuse silent overwrite."""
    if path.exists() and not overwrite:
        raise FileExistsError(f"Output already exists: {path}. Use --overwrite explicitly.")
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            "w",
            encoding="utf-8",
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            temporary_path = Path(handle.name)
            json.dump(payload, handle, indent=2)
            handle.write("\n")
        os.replace(temporary_path, path)
        temporary_path = None
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)
