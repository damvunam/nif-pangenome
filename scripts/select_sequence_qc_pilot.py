"""Select a deterministic 10-genome pilot for sequence-based QC."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import tempfile
from collections.abc import Callable, Mapping, Sequence
from pathlib import Path
from typing import Any


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT_PATH = (
    REPOSITORY_ROOT / "data/metadata/bradyrhizobium_metadata_qc.tsv"
)
DEFAULT_POLICY_PATH = REPOSITORY_ROOT / "config/metadata_qc_policy.json"
DEFAULT_OUTPUT_PATH = (
    REPOSITORY_ROOT
    / "data/metadata/bradyrhizobium_sequence_qc_pilot_10.tsv"
)

REQUIRED_INPUT_FIELDS = {
    "canonical_id",
    "pair_status",
    "selected_accession",
    "selected_source",
    "selected_status",
    "selection_reason",
    "organism_name",
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
}

PILOT_FIELDS = [
    "pilot_order",
    "pilot_role",
    "selection_rule",
    "canonical_id",
    "selected_accession",
    "selected_source",
    "selected_status",
    "selection_reason",
    "pair_status",
    "organism_name",
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
    "source_metadata_qc_sha256",
    "source_policy_sha256",
]


class PilotSelectionError(ValueError):
    """Raised when the metadata cannot satisfy the locked pilot design."""


def sha256_file(path: Path) -> str:
    """Return the SHA-256 digest of a file."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_metadata_qc(path: Path) -> list[dict[str, str]]:
    """Load metadata-QC rows and enforce the pilot input contract."""
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        fields = set(reader.fieldnames or [])
    missing = sorted(REQUIRED_INPUT_FIELDS - fields)
    if missing:
        raise PilotSelectionError(f"Metadata QC is missing columns: {missing!r}.")
    if not rows:
        raise PilotSelectionError("Metadata QC input is empty.")
    canonical_ids = [row["canonical_id"] for row in rows]
    accessions = [row["selected_accession"] for row in rows]
    if len(canonical_ids) != len(set(canonical_ids)):
        raise PilotSelectionError("Metadata QC contains duplicate canonical IDs.")
    if len(accessions) != len(set(accessions)):
        raise PilotSelectionError("Metadata QC contains duplicate selected accessions.")
    if any(row["selected_status"] != "current" for row in rows):
        raise PilotSelectionError("All pilot candidates must be current assemblies.")
    return rows


def load_contiguity_policy(path: Path) -> dict[str, float]:
    """Load the locked metadata contiguity thresholds used by the pilot."""
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
        contiguity = payload["metadata_contiguity"]
        values = {
            "synteny_max_contigs": float(contiguity["synteny_max_contigs"]),
            "synteny_min_contig_n50_bp": float(
                contiguity["synteny_min_contig_n50_bp"]
            ),
        }
    except (OSError, json.JSONDecodeError, KeyError, TypeError, ValueError) as exc:
        raise PilotSelectionError(f"Invalid metadata-QC policy {path}: {exc}.") from exc
    if any(value <= 0 for value in values.values()):
        raise PilotSelectionError("Synteny contiguity thresholds must be positive.")
    return values


def _flags(row: Mapping[str, str]) -> frozenset[str]:
    return frozenset(filter(None, row["flags"].split(";")))


def _float(row: Mapping[str, str], field: str) -> float:
    try:
        return float(row[field])
    except (TypeError, ValueError) as exc:
        raise PilotSelectionError(
            f"{row['selected_accession']} has invalid {field}: {row[field]!r}."
        ) from exc


def _int(row: Mapping[str, str], field: str) -> int:
    try:
        return int(row[field])
    except (TypeError, ValueError) as exc:
        raise PilotSelectionError(
            f"{row['selected_accession']} has invalid {field}: {row[field]!r}."
        ) from exc


def _center_score(row: Mapping[str, str]) -> tuple[float, str]:
    return (
        abs(_float(row, "genome_size_modified_z"))
        + abs(_float(row, "gc_percent_modified_z")),
        row["selected_accession"],
    )


def _choose(
    rows: Sequence[dict[str, str]],
    *,
    role: str,
    rule: str,
    predicate: Callable[[dict[str, str]], bool],
    sort_key: Callable[[dict[str, str]], Any],
    used: set[str],
) -> tuple[str, str, dict[str, str]]:
    candidates = [
        row
        for row in rows
        if row["selected_accession"] not in used and predicate(row)
    ]
    if not candidates:
        raise PilotSelectionError(f"No eligible candidate for pilot role {role!r}.")
    selected = min(candidates, key=sort_key)
    used.add(selected["selected_accession"])
    return role, rule, selected


def select_pilot_rows(
    rows: Sequence[dict[str, str]], contiguity: Mapping[str, float]
) -> list[tuple[str, str, dict[str, str]]]:
    """Select ten non-overlapping genomes using the approved pilot strata."""
    used: set[str] = set()
    selected: list[tuple[str, str, dict[str, str]]] = []

    clean_refseq = lambda row: (
        row["selected_source"] == "REFSEQ"
        and not _flags(row)
        and row["metadata_decision"] == "retain_candidate"
        and row["pangenome_contiguity_class"] == "primary_contiguity_candidate"
        and row["synteny_class"] == "synteny_candidate"
    )
    for index in (1, 2):
        selected.append(
            _choose(
                rows,
                role=f"normal_refseq_near_median_{index}",
                rule="Clean RefSeq and minimum combined absolute metadata z-score",
                predicate=clean_refseq,
                sort_key=_center_score,
                used=used,
            )
        )

    selected.append(
        _choose(
            rows,
            role="normal_genbank_singleton_near_median",
            rule="Clean GenBank singleton and minimum combined absolute metadata z-score",
            predicate=lambda row: (
                row["selection_reason"] == "current_genbank_singleton"
                and not _flags(row)
                and row["metadata_decision"] == "retain_candidate"
                and row["pangenome_contiguity_class"]
                == "primary_contiguity_candidate"
                and row["synteny_class"] == "synteny_candidate"
            ),
            sort_key=_center_score,
            used=used,
        )
    )

    selected.append(
        _choose(
            rows,
            role="fragmented_extreme",
            rule="Isolated pangenome/synteny fragmentation and maximum contig count",
            predicate=lambda row: _flags(row)
            == frozenset({"pangenome_fragmented", "synteny_fragmented"}),
            sort_key=lambda row: (
                -_int(row, "number_of_contigs"),
                _int(row, "contig_n50"),
                row["selected_accession"],
            ),
            used=used,
        )
    )

    selected.append(
        _choose(
            rows,
            role="genome_size_outlier_isolated",
            rule="Genome-size-only outlier with maximum absolute modified z-score",
            predicate=lambda row: _flags(row) == frozenset({"genome_size_outlier"}),
            sort_key=lambda row: (
                -abs(_float(row, "genome_size_modified_z")),
                row["selected_accession"],
            ),
            used=used,
        )
    )

    selected.append(
        _choose(
            rows,
            role="gc_percent_outlier_isolated",
            rule="GC-only outlier with maximum absolute modified z-score",
            predicate=lambda row: _flags(row) == frozenset({"gc_percent_outlier"}),
            sort_key=lambda row: (
                -abs(_float(row, "gc_percent_modified_z")),
                row["selected_accession"],
            ),
            used=used,
        )
    )

    selected.append(
        _choose(
            rows,
            role="missing_gc_percent_isolated",
            rule="Missing-GC-only case nearest the genome-size median",
            predicate=lambda row: _flags(row) == frozenset({"missing_gc_percent"}),
            sort_key=lambda row: (
                abs(_float(row, "genome_size_modified_z")),
                row["selected_accession"],
            ),
            used=used,
        )
    )

    selected.append(
        _choose(
            rows,
            role="refseq_to_genbank_fallback_isolated",
            rule="Fallback-only case nearest both metadata medians",
            predicate=lambda row: _flags(row)
            == frozenset({"suppressed_or_unavailable_refseq_fallback"}),
            sort_key=_center_score,
            used=used,
        )
    )

    selected.append(
        _choose(
            rows,
            role="broken_pair",
            rule="Broken-pair case nearest both metadata medians",
            predicate=lambda row: "broken_pair" in _flags(row),
            sort_key=_center_score,
            used=used,
        )
    )

    max_contigs = contiguity["synteny_max_contigs"]
    min_n50 = contiguity["synteny_min_contig_n50_bp"]

    def synteny_boundary_distance(row: Mapping[str, str]) -> tuple[float, str]:
        contig_excess = max(0.0, _int(row, "number_of_contigs") - max_contigs)
        n50_shortfall = max(0.0, min_n50 - _int(row, "contig_n50"))
        return (
            contig_excess / max_contigs + n50_shortfall / min_n50,
            row["selected_accession"],
        )

    selected.append(
        _choose(
            rows,
            role="synteny_threshold_edge",
            rule="Clean primary candidate with minimum normalized synteny-threshold miss",
            predicate=lambda row: (
                _flags(row) == frozenset({"synteny_fragmented"})
                and row["pangenome_contiguity_class"]
                == "primary_contiguity_candidate"
                and row["synteny_class"] == "presence_absence_only"
            ),
            sort_key=synteny_boundary_distance,
            used=used,
        )
    )

    if len(selected) != 10 or len(used) != 10:
        raise PilotSelectionError("Pilot selection must contain 10 unique genomes.")
    return selected


def build_output_rows(
    selected: Sequence[tuple[str, str, dict[str, str]]],
    *,
    metadata_qc_sha256: str,
    policy_sha256: str,
) -> list[dict[str, str | int]]:
    """Build the declared pilot TSV rows with input provenance."""
    output: list[dict[str, str | int]] = []
    copied_fields = PILOT_FIELDS[3:-2]
    for order, (role, rule, source) in enumerate(selected, start=1):
        row: dict[str, str | int] = {
            "pilot_order": order,
            "pilot_role": role,
            "selection_rule": rule,
        }
        row.update({field: source[field] for field in copied_fields})
        row["source_metadata_qc_sha256"] = metadata_qc_sha256
        row["source_policy_sha256"] = policy_sha256
        output.append(row)
    return output


def write_pilot_tsv(
    rows: Sequence[Mapping[str, Any]], path: Path, *, overwrite: bool
) -> None:
    """Write the pilot TSV atomically and refuse silent overwrite."""
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
                fieldnames=PILOT_FIELDS,
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
    parser.add_argument("--input-path", type=Path, default=DEFAULT_INPUT_PATH)
    parser.add_argument("--policy-path", type=Path, default=DEFAULT_POLICY_PATH)
    parser.add_argument("--output-path", type=Path, default=DEFAULT_OUTPUT_PATH)
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Explicitly replace an existing pilot TSV atomically.",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    """Select the pilot, write its manifest, and report that no QC was run."""
    args = parse_args(argv)
    rows = load_metadata_qc(args.input_path)
    contiguity = load_contiguity_policy(args.policy_path)
    selected = select_pilot_rows(rows, contiguity)
    output_rows = build_output_rows(
        selected,
        metadata_qc_sha256=sha256_file(args.input_path),
        policy_sha256=sha256_file(args.policy_path),
    )
    write_pilot_tsv(output_rows, args.output_path, overwrite=args.overwrite)
    print(f"Pilot genomes selected: {len(output_rows)}")
    for row in output_rows:
        print(f"{row['pilot_order']:02d} {row['selected_accession']} {row['pilot_role']}")
    print("Sequence QC executed: NO")
    print(f"Pilot TSV: {args.output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
