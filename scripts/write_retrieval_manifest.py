"""Recompute raw-file integrity and empirical counts in a retrieval manifest.

Retrieval provenance is preserved from an existing manifest because timestamps,
tool versions, and operator identity cannot be reconstructed from raw content.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import tempfile
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_RAW_PATH = REPOSITORY_ROOT / "data/raw/bradyrhizobium_genome_summary.jsonl"
DEFAULT_MANIFEST_PATH = REPOSITORY_ROOT / "data/metadata/retrieval_manifest.json"


class ManifestError(ValueError):
    """Raised when a manifest or raw snapshot violates the data contract."""


def sha256sum(path: Path) -> str:
    """Return the SHA-256 digest of a file."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_records(raw_path: Path) -> list[dict[str, Any]]:
    """Read JSON Lines records and reject duplicates or malformed accessions."""
    records: list[dict[str, Any]] = []
    seen: set[str] = set()
    with raw_path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                raise ManifestError(f"Blank JSONL record at line {line_number}.")
            try:
                record = json.loads(line)
            except json.JSONDecodeError as exc:
                raise ManifestError(
                    f"Invalid JSON at line {line_number}: {exc.msg}."
                ) from exc
            if not isinstance(record, dict):
                raise ManifestError(f"Record at line {line_number} is not an object.")
            accession = record.get("accession")
            if not isinstance(accession, str) or not accession:
                raise ManifestError(f"Missing accession at line {line_number}.")
            if accession in seen:
                raise ManifestError(f"Duplicate accession {accession!r}.")
            seen.add(accession)
            records.append(record)
    if not records:
        raise ManifestError(f"No records found in {raw_path}.")
    return records


def build_manifest(
    raw_path: Path,
    records: Sequence[Mapping[str, Any]],
    retrieval: Mapping[str, Any],
) -> dict[str, Any]:
    """Build content-derived manifest fields while preserving retrieval facts."""
    source_counts: dict[str, int] = {}
    for record in records:
        source = str(record.get("source_database", "MISSING"))
        source_counts[source] = source_counts.get(source, 0) + 1
    paired_count = sum(bool(record.get("paired_accession")) for record in records)
    return {
        "retrieval": dict(retrieval),
        "raw_file": {
            "path": "data/raw/bradyrhizobium_genome_summary.jsonl",
            "size_bytes": raw_path.stat().st_size,
            "sha256": sha256sum(raw_path),
            "git_tracked": False,
            "note": "Excluded from git via .gitignore; integrity tracked here via checksum.",
        },
        "empirical_counts": {
            "total_primary_records": len(records),
            "distinct_accessions": len({record["accession"] for record in records}),
            "source_database_counts": source_counts,
            "records_with_paired_accession": paired_count,
            "records_without_paired_accession": len(records) - paired_count,
        },
    }


def load_retrieval_provenance(path: Path) -> Mapping[str, Any]:
    """Load the non-reconstructable retrieval section from a prior manifest."""
    try:
        manifest = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ManifestError(f"Cannot read provenance manifest {path}: {exc}.") from exc
    retrieval = manifest.get("retrieval") if isinstance(manifest, dict) else None
    if not isinstance(retrieval, dict) or not retrieval:
        raise ManifestError(f"Manifest {path} has no retrieval provenance section.")
    return retrieval


def write_manifest(manifest: Mapping[str, Any], path: Path, *, overwrite: bool) -> None:
    """Write the manifest atomically and refuse silent replacement."""
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
            json.dump(manifest, handle, indent=2)
            handle.write("\n")
        os.replace(temporary_path, path)
        temporary_path = None
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-path", type=Path, default=DEFAULT_RAW_PATH)
    parser.add_argument("--output-path", type=Path, default=DEFAULT_MANIFEST_PATH)
    parser.add_argument(
        "--provenance-from",
        type=Path,
        default=DEFAULT_MANIFEST_PATH,
        help="Existing manifest whose retrieval section will be preserved.",
    )
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    """Recompute and safely write the retrieval manifest."""
    args = parse_args(argv)
    retrieval = load_retrieval_provenance(args.provenance_from)
    records = load_records(args.raw_path)
    manifest = build_manifest(args.raw_path, records, retrieval)
    write_manifest(manifest, args.output_path, overwrite=args.overwrite)
    print(f"Manifest written to {args.output_path}")
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
