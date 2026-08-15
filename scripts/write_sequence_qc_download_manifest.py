"""Validate the pilot genome package and write its download manifest."""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import re
import tempfile
from collections.abc import Mapping, Sequence
from datetime import datetime
from pathlib import Path
from typing import Any


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_PILOT_PATH = (
    REPOSITORY_ROOT
    / "data/metadata/bradyrhizobium_sequence_qc_pilot_10.tsv"
)
DEFAULT_RAW_ROOT = REPOSITORY_ROOT / "data/raw/sequence_qc_pilot_10"
DEFAULT_OUTPUT_PATH = (
    REPOSITORY_ROOT / "data/metadata/bradyrhizobium_download_manifest.tsv"
)

MANIFEST_SCHEMA_VERSION = 1
EXPECTED_PILOT_COUNT = 10
PILOT_FIELDS = {
    "pilot_order",
    "pilot_role",
    "canonical_id",
    "selected_accession",
}
MANIFEST_FIELDS = [
    "manifest_schema_version",
    "pilot_order",
    "pilot_role",
    "canonical_id",
    "selected_accession",
    "fasta_package_path",
    "fasta_size_bytes",
    "fasta_sha256",
    "ncbi_package_size_bytes",
    "ncbi_package_sha256",
    "ncbi_md5_validation",
    "datasets_version",
    "retrieval_timestamp_local",
    "source_git_commit",
    "source_pilot_sha256",
]
SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")
GIT_COMMIT_PATTERN = re.compile(r"^[0-9a-f]{40}$")


class DownloadManifestError(ValueError):
    """Raised when downloaded pilot data violate the provenance contract."""


def file_digest(path: Path, algorithm: str = "sha256") -> str:
    """Return a hexadecimal digest for a file."""
    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_pilot(path: Path) -> list[dict[str, str]]:
    """Load the locked ten-genome pilot in declared order."""
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        fields = set(reader.fieldnames or [])
    missing = sorted(PILOT_FIELDS - fields)
    if missing:
        raise DownloadManifestError(f"Pilot TSV is missing columns: {missing!r}.")
    if len(rows) != EXPECTED_PILOT_COUNT:
        raise DownloadManifestError(
            f"Pilot TSV must contain {EXPECTED_PILOT_COUNT} rows, found {len(rows)}."
        )
    expected_orders = [str(index) for index in range(1, EXPECTED_PILOT_COUNT + 1)]
    if [row["pilot_order"] for row in rows] != expected_orders:
        raise DownloadManifestError("Pilot orders must be consecutive from 1 through 10.")
    accessions = [row["selected_accession"] for row in rows]
    canonical_ids = [row["canonical_id"] for row in rows]
    if len(accessions) != len(set(accessions)):
        raise DownloadManifestError("Pilot TSV contains duplicate accessions.")
    if len(canonical_ids) != len(set(canonical_ids)):
        raise DownloadManifestError("Pilot TSV contains duplicate canonical IDs.")
    return rows


def _read_single_line(path: Path, label: str) -> str:
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError as exc:
        raise DownloadManifestError(f"Cannot read {label} file {path}: {exc}.") from exc
    if len(lines) != 1 or not lines[0].strip():
        raise DownloadManifestError(f"{label} file must contain exactly one value.")
    return lines[0].strip()


def _parse_checksum_file(path: Path, label: str) -> tuple[str, str]:
    line = _read_single_line(path, label)
    parts = line.split(maxsplit=1)
    if len(parts) != 2 or not SHA256_PATTERN.fullmatch(parts[0]):
        raise DownloadManifestError(f"Invalid SHA-256 record in {path}.")
    return parts[0], parts[1].lstrip("* ")


def _load_accessions(path: Path) -> list[str]:
    try:
        values = path.read_text(encoding="utf-8").splitlines()
    except OSError as exc:
        raise DownloadManifestError(f"Cannot read accession file {path}: {exc}.") from exc
    if any(not value for value in values):
        raise DownloadManifestError("Downloaded accession list contains blank values.")
    if len(values) != len(set(values)):
        raise DownloadManifestError("Downloaded accession list contains duplicates.")
    return values


def _load_fasta_checksums(path: Path) -> dict[str, str]:
    checksums: dict[str, str] = {}
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError as exc:
        raise DownloadManifestError(f"Cannot read FASTA checksum file {path}: {exc}.") from exc
    for line_number, line in enumerate(lines, start=1):
        parts = line.split(maxsplit=1)
        if len(parts) != 2 or not SHA256_PATTERN.fullmatch(parts[0]):
            raise DownloadManifestError(
                f"Invalid FASTA SHA-256 record at line {line_number}."
            )
        relative_path = parts[1].lstrip("* ")
        if relative_path in checksums:
            raise DownloadManifestError(
                f"Duplicate FASTA checksum path {relative_path!r}."
            )
        checksums[relative_path] = parts[0]
    if len(checksums) != EXPECTED_PILOT_COUNT:
        raise DownloadManifestError(
            f"Expected {EXPECTED_PILOT_COUNT} FASTA checksums, found {len(checksums)}."
        )
    return checksums


def _validate_ncbi_md5(unpacked_root: Path) -> None:
    checksum_path = unpacked_root / "md5sum.txt"
    try:
        lines = checksum_path.read_text(encoding="utf-8").splitlines()
    except OSError as exc:
        raise DownloadManifestError(f"Cannot read NCBI MD5 file: {exc}.") from exc
    if not lines:
        raise DownloadManifestError("NCBI MD5 file is empty.")
    for line_number, line in enumerate(lines, start=1):
        parts = line.split(maxsplit=1)
        if len(parts) != 2 or not re.fullmatch(r"[0-9a-fA-F]{32}", parts[0]):
            raise DownloadManifestError(
                f"Invalid NCBI MD5 record at line {line_number}."
            )
        relative_path = parts[1].lstrip("* ")
        target = unpacked_root / relative_path
        if not target.is_file():
            raise DownloadManifestError(f"NCBI package file is missing: {relative_path}.")
        if file_digest(target, "md5").lower() != parts[0].lower():
            raise DownloadManifestError(
                f"NCBI MD5 mismatch for package file {relative_path}."
            )


def _load_provenance(raw_root: Path, pilot_path: Path) -> dict[str, str | int]:
    source_pilot_sha256, _ = _parse_checksum_file(
        raw_root / "source_pilot.sha256", "source pilot checksum"
    )
    actual_pilot_sha256 = file_digest(pilot_path)
    if source_pilot_sha256 != actual_pilot_sha256:
        raise DownloadManifestError("Pilot TSV SHA-256 does not match download provenance.")

    package_path = raw_root / "ncbi_dataset.zip"
    package_sha256, package_name = _parse_checksum_file(
        raw_root / "ncbi_dataset.zip.sha256", "NCBI package checksum"
    )
    if package_name != package_path.name:
        raise DownloadManifestError("NCBI package checksum names an unexpected file.")
    if not package_path.is_file():
        raise DownloadManifestError(f"NCBI package is missing: {package_path}.")
    if file_digest(package_path) != package_sha256:
        raise DownloadManifestError("NCBI package SHA-256 mismatch.")

    version_line = _read_single_line(
        raw_root / "datasets_version.txt", "NCBI Datasets version"
    )
    version_match = re.fullmatch(r"datasets version:\s*(\S+)", version_line)
    if version_match is None:
        raise DownloadManifestError("Invalid NCBI Datasets version record.")

    timestamp = _read_single_line(
        raw_root / "retrieval_timestamp_local.txt", "retrieval timestamp"
    )
    try:
        datetime.fromisoformat(timestamp)
    except ValueError as exc:
        raise DownloadManifestError("Retrieval timestamp is not valid ISO 8601.") from exc

    source_git_commit = _read_single_line(
        raw_root / "source_git_commit.txt", "source Git commit"
    )
    if not GIT_COMMIT_PATTERN.fullmatch(source_git_commit):
        raise DownloadManifestError("Source Git commit must be a full 40-character SHA.")

    return {
        "ncbi_package_size_bytes": package_path.stat().st_size,
        "ncbi_package_sha256": package_sha256,
        "datasets_version": version_match.group(1),
        "retrieval_timestamp_local": timestamp,
        "source_git_commit": source_git_commit,
        "source_pilot_sha256": source_pilot_sha256,
    }


def build_download_manifest(
    pilot_rows: Sequence[Mapping[str, str]],
    pilot_path: Path,
    raw_root: Path,
) -> list[dict[str, str | int]]:
    """Validate the downloaded package and build deterministic manifest rows."""
    expected_accessions = [row["selected_accession"] for row in pilot_rows]
    downloaded_accessions = _load_accessions(raw_root / "accessions.txt")
    if downloaded_accessions != expected_accessions:
        raise DownloadManifestError(
            "Downloaded accession order does not match the locked pilot TSV."
        )

    provenance = _load_provenance(raw_root, pilot_path)
    unpacked_root = raw_root / "unpacked"
    data_root = unpacked_root / "ncbi_dataset/data"
    if not data_root.is_dir():
        raise DownloadManifestError(f"NCBI data directory is missing: {data_root}.")
    _validate_ncbi_md5(unpacked_root)
    recorded_fasta_checksums = _load_fasta_checksums(
        raw_root / "genome_fasta_sha256.txt"
    )

    downloaded_directories = sorted(
        path.name for path in data_root.iterdir() if path.is_dir()
    )
    if downloaded_directories != sorted(expected_accessions):
        raise DownloadManifestError(
            "Downloaded assembly directories do not match the locked pilot accessions."
        )

    output: list[dict[str, str | int]] = []
    observed_checksum_paths: set[str] = set()
    for pilot_row in pilot_rows:
        accession = pilot_row["selected_accession"]
        accession_root = data_root / accession
        fasta_paths = sorted(accession_root.glob("*_genomic.fna"))
        if len(fasta_paths) != 1:
            raise DownloadManifestError(
                f"Expected one genomic FASTA for {accession}, found {len(fasta_paths)}."
            )
        fasta_path = fasta_paths[0]
        if fasta_path.stat().st_size == 0:
            raise DownloadManifestError(f"Genomic FASTA is empty for {accession}.")
        with fasta_path.open("rb") as handle:
            if handle.read(1) != b">":
                raise DownloadManifestError(
                    f"Genomic FASTA has an invalid first byte for {accession}."
                )

        package_path = fasta_path.relative_to(unpacked_root).as_posix()
        fasta_sha256 = file_digest(fasta_path)
        recorded_sha256 = recorded_fasta_checksums.get(package_path)
        if recorded_sha256 != fasta_sha256:
            raise DownloadManifestError(
                f"Recorded FASTA SHA-256 mismatch for {accession}."
            )
        observed_checksum_paths.add(package_path)

        manifest_row: dict[str, str | int] = {
            "manifest_schema_version": MANIFEST_SCHEMA_VERSION,
            "pilot_order": pilot_row["pilot_order"],
            "pilot_role": pilot_row["pilot_role"],
            "canonical_id": pilot_row["canonical_id"],
            "selected_accession": accession,
            "fasta_package_path": package_path,
            "fasta_size_bytes": fasta_path.stat().st_size,
            "fasta_sha256": fasta_sha256,
            "ncbi_md5_validation": "pass",
        }
        manifest_row.update(provenance)
        output.append(manifest_row)

    if observed_checksum_paths != set(recorded_fasta_checksums):
        raise DownloadManifestError("FASTA checksum paths do not match selected genomes.")
    return output


def write_manifest_tsv(
    rows: Sequence[Mapping[str, Any]], path: Path, *, overwrite: bool
) -> None:
    """Write the manifest atomically and refuse silent overwrite."""
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
                fieldnames=MANIFEST_FIELDS,
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


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pilot-path", type=Path, default=DEFAULT_PILOT_PATH)
    parser.add_argument("--raw-root", type=Path, default=DEFAULT_RAW_ROOT)
    parser.add_argument("--output-path", type=Path, default=DEFAULT_OUTPUT_PATH)
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Explicitly replace an existing manifest atomically.",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    """Validate the package and write the tracked download manifest."""
    args = parse_args(argv)
    pilot_rows = load_pilot(args.pilot_path)
    manifest_rows = build_download_manifest(
        pilot_rows, args.pilot_path, args.raw_root
    )
    write_manifest_tsv(manifest_rows, args.output_path, overwrite=args.overwrite)
    print(f"Download manifest rows: {len(manifest_rows)}")
    print("NCBI package validation: PASS")
    print("FASTA SHA-256 validation: PASS")
    print("Sequence QC executed: NO")
    print(f"Manifest TSV: {args.output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
