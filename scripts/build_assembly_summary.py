"""Build a deterministic canonical Bradyrhizobium assembly summary.

GCA/GCF pairs are collapsed to one row. A present RefSeq record is used as
the representative; otherwise, the available singleton record is retained.
No biological quality threshold is applied by this script.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import tempfile
from collections import Counter
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_RAW_PATH = REPOSITORY_ROOT / "data/raw/bradyrhizobium_genome_summary.jsonl"
DEFAULT_OUTPUT_PATH = REPOSITORY_ROOT / "data/metadata/bradyrhizobium_assembly_summary.tsv"

FIELDS = [
    "canonical_id",
    "representative_accession",
    "representative_source",
    "gcf_accession",
    "gca_accession",
    "pair_status",
    "missing_partner_accession",
    "organism_name",
    "tax_id",
    "strain",
    "assembly_level",
    "assembly_status",
    "refseq_category",
    "genome_size_bp",
    "gc_percent",
    "number_of_contigs",
    "number_of_scaffolds",
    "contig_n50",
    "scaffold_n50",
    "bioproject_accession",
    "biosample_accession",
    "submitter",
    "release_date",
]

SOURCE_REFSEQ = "SOURCE_DATABASE_REFSEQ"
SOURCE_GENBANK = "SOURCE_DATABASE_GENBANK"


class AssemblySummaryError(ValueError):
    """Raised when raw assembly metadata violates the summary data contract."""


def source_label(record: Mapping[str, Any]) -> str:
    """Return REFSEQ or GENBANK while validating the accession prefix."""
    accession = record.get("accession")
    source = record.get("source_database")
    if not isinstance(accession, str) or not accession:
        raise AssemblySummaryError("Every record must have a non-empty accession.")
    if source == SOURCE_REFSEQ and accession.startswith("GCF_"):
        return "REFSEQ"
    if source == SOURCE_GENBANK and accession.startswith("GCA_"):
        return "GENBANK"
    raise AssemblySummaryError(
        f"Accession/source mismatch for {accession!r}: {source!r}."
    )


def load_records(raw_path: Path) -> dict[str, dict[str, Any]]:
    """Load JSON Lines records and reject blank, malformed, or duplicate entries."""
    records: dict[str, dict[str, Any]] = {}
    with raw_path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                raise AssemblySummaryError(
                    f"Blank JSONL record at {raw_path}:{line_number}."
                )
            try:
                record = json.loads(line)
            except json.JSONDecodeError as exc:
                raise AssemblySummaryError(
                    f"Invalid JSON at {raw_path}:{line_number}: {exc.msg}."
                ) from exc
            if not isinstance(record, dict):
                raise AssemblySummaryError(
                    f"Record at {raw_path}:{line_number} is not a JSON object."
                )
            accession = record.get("accession")
            source_label(record)
            if accession in records:
                raise AssemblySummaryError(
                    f"Duplicate accession {accession!r} at line {line_number}."
                )
            records[accession] = record
    if not records:
        raise AssemblySummaryError(f"No records found in {raw_path}.")
    return records


def extract_row(
    record: Mapping[str, Any],
    *,
    gcf_accession: str,
    gca_accession: str,
    pair_status: str,
    missing_partner_accession: str = "",
) -> dict[str, Any]:
    """Extract the stable TSV schema from one representative record."""
    organism = record.get("organism", {}) or {}
    assembly_info = record.get("assembly_info", {}) or {}
    assembly_stats = record.get("assembly_stats", {}) or {}
    biosample = assembly_info.get("biosample", {}) or {}

    return {
        "canonical_id": "",
        "representative_accession": record.get("accession"),
        "representative_source": source_label(record),
        "gcf_accession": gcf_accession,
        "gca_accession": gca_accession,
        "pair_status": pair_status,
        "missing_partner_accession": missing_partner_accession,
        "organism_name": organism.get("organism_name", ""),
        "tax_id": organism.get("tax_id", ""),
        "strain": (organism.get("infraspecific_names", {}) or {}).get("strain", ""),
        "assembly_level": assembly_info.get("assembly_level", ""),
        "assembly_status": assembly_info.get("assembly_status", ""),
        "refseq_category": assembly_info.get("refseq_category", ""),
        "genome_size_bp": assembly_stats.get("total_sequence_length", ""),
        "gc_percent": assembly_stats.get("gc_percent", ""),
        "number_of_contigs": assembly_stats.get("number_of_contigs", ""),
        "number_of_scaffolds": assembly_stats.get("number_of_scaffolds", ""),
        "contig_n50": assembly_stats.get("contig_n50", ""),
        "scaffold_n50": assembly_stats.get("scaffold_n50", ""),
        "bioproject_accession": assembly_info.get("bioproject_accession", ""),
        "biosample_accession": biosample.get("accession", ""),
        "submitter": assembly_info.get("submitter", ""),
        "release_date": assembly_info.get("release_date", ""),
    }


def build_rows(records: Mapping[str, Mapping[str, Any]]) -> list[dict[str, Any]]:
    """Collapse records into validated canonical groups in stable accession order."""
    seen: set[str] = set()
    rows: list[dict[str, Any]] = []

    for accession in sorted(records):
        if accession in seen:
            continue
        record = records[accession]
        partner_accession = record.get("paired_accession")
        if partner_accession is not None and not isinstance(partner_accession, str):
            raise AssemblySummaryError(
                f"paired_accession for {accession!r} must be a string or null."
            )
        partner = records.get(partner_accession) if partner_accession else None

        if partner is not None:
            if partner.get("paired_accession") != accession:
                raise AssemblySummaryError(
                    f"Non-reciprocal pair: {accession!r} points to "
                    f"{partner_accession!r}, but the partner points to "
                    f"{partner.get('paired_accession')!r}."
                )
            labels = {source_label(record), source_label(partner)}
            if labels != {"REFSEQ", "GENBANK"}:
                raise AssemblySummaryError(
                    f"Pair {accession!r}/{partner_accession!r} does not contain "
                    "exactly one RefSeq and one GenBank record."
                )
            refseq_record = record if source_label(record) == "REFSEQ" else partner
            genbank_record = record if source_label(record) == "GENBANK" else partner
            rows.append(
                extract_row(
                    refseq_record,
                    gcf_accession=str(refseq_record["accession"]),
                    gca_accession=str(genbank_record["accession"]),
                    pair_status="paired",
                )
            )
            seen.update({accession, str(partner_accession)})
            continue

        label = source_label(record)
        if partner_accession:
            expected_prefix = "GCA_" if label == "REFSEQ" else "GCF_"
            if not partner_accession.startswith(expected_prefix):
                raise AssemblySummaryError(
                    f"Broken-pair partner for {accession!r} must start with "
                    f"{expected_prefix!r}, not {partner_accession!r}."
                )
        rows.append(
            extract_row(
                record,
                gcf_accession=accession if label == "REFSEQ" else "",
                gca_accession=accession if label == "GENBANK" else "",
                pair_status="broken_pair" if partner_accession else "unpaired",
                missing_partner_accession=str(partner_accession or ""),
            )
        )
        seen.add(accession)

    if seen != set(records):
        missing = sorted(set(records) - seen)
        raise AssemblySummaryError(
            f"Not all input records were consumed: {missing[:5]!r}."
        )

    rows.sort(key=lambda row: (row["gca_accession"] or row["gcf_accession"]))
    for index, row in enumerate(rows, start=1):
        row["canonical_id"] = f"CANON_{index:05d}"
    return rows


def write_rows(
    rows: Sequence[Mapping[str, Any]], output_path: Path, *, overwrite: bool
) -> None:
    """Write TSV atomically and refuse silent replacement by default."""
    if output_path.exists() and not overwrite:
        raise FileExistsError(
            f"Output already exists: {output_path}. Use --overwrite explicitly."
        )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            "w",
            encoding="utf-8",
            newline="",
            dir=output_path.parent,
            prefix=f".{output_path.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            temporary_path = Path(handle.name)
            writer = csv.DictWriter(
                handle,
                fieldnames=FIELDS,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temporary_path, output_path)
        temporary_path = None
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-path", type=Path, default=DEFAULT_RAW_PATH)
    parser.add_argument("--output-path", type=Path, default=DEFAULT_OUTPUT_PATH)
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Explicitly replace an existing output file atomically.",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    """Build the assembly summary and print a concise audit breakdown."""
    args = parse_args(argv)
    records = load_records(args.raw_path)
    rows = build_rows(records)
    write_rows(rows, args.output_path, overwrite=args.overwrite)
    print(f"Wrote {len(rows)} canonical assembly rows to {args.output_path}")
    print(f"Total input records consumed: {len(records)}")
    print("pair_status breakdown:", Counter(row["pair_status"] for row in rows))
    print(
        "representative_source breakdown:",
        Counter(row["representative_source"] for row in rows),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())