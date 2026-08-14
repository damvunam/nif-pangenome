"""Reusable metadata-only quality control for canonical genome assemblies."""

from __future__ import annotations

import math
import statistics
from collections import Counter
from collections.abc import Mapping, Sequence
from typing import Any

try:
    from .qc_io import (
        AUDIT_FIELDS,
        MetadataQCError,
        load_policy,
        load_raw_records,
        load_summary,
        sha256_file,
        write_json_atomic,
        write_tsv_atomic,
    )
except ImportError:  # Direct execution from scripts/.
    from qc_io import (
        AUDIT_FIELDS,
        MetadataQCError,
        load_policy,
        load_raw_records,
        load_summary,
        sha256_file,
        write_json_atomic,
        write_tsv_atomic,
    )


def record_status(record: Mapping[str, Any]) -> str:
    """Return a normalized NCBI assembly status."""
    assembly_info = record.get("assembly_info", {}) or {}
    return str(assembly_info.get("assembly_status", "")).strip().lower()


def record_source(record: Mapping[str, Any]) -> str:
    """Return REFSEQ or GENBANK while validating the accession/source pair."""
    accession = str(record.get("accession", ""))
    source = record.get("source_database")
    if accession.startswith("GCF_") and source == "SOURCE_DATABASE_REFSEQ":
        return "REFSEQ"
    if accession.startswith("GCA_") and source == "SOURCE_DATABASE_GENBANK":
        return "GENBANK"
    raise MetadataQCError(
        f"Accession/source mismatch for {accession!r}: {source!r}."
    )


def validate_member_coverage(
    summary_rows: Sequence[Mapping[str, str]],
    raw_records: Mapping[str, Mapping[str, Any]],
) -> None:
    """Require every accession represented in the TSV to exist in raw metadata."""
    represented = {
        accession
        for row in summary_rows
        for accession in (row["gcf_accession"], row["gca_accession"])
        if accession
    }
    missing = sorted(represented - set(raw_records))
    if missing:
        raise MetadataQCError(
            f"Raw metadata is missing represented accessions: {missing[:5]!r}."
        )
    extra = sorted(set(raw_records) - represented)
    if extra:
        raise MetadataQCError(
            f"Raw metadata contains accessions absent from the canonical summary: "
            f"{extra[:5]!r}."
        )


def select_representative(
    row: Mapping[str, str],
    raw_records: Mapping[str, Mapping[str, Any]],
    current_statuses: set[str],
) -> tuple[Mapping[str, Any] | None, str]:
    """Select a current representative and return a traceable reason."""
    members = [
        raw_records[accession]
        for accession in (row["gcf_accession"], row["gca_accession"])
        if accession
    ]
    current = [record for record in members if record_status(record) in current_statuses]
    current_refseq = [record for record in current if record_source(record) == "REFSEQ"]
    if current_refseq:
        return current_refseq[0], "current_refseq_preferred"
    current_genbank = [record for record in current if record_source(record) == "GENBANK"]
    if current_genbank:
        selected = current_genbank[0]
        original = row["representative_accession"]
        if row["pair_status"] == "broken_pair":
            reason = "broken_pair_current_genbank"
        elif original != selected["accession"]:
            reason = "fallback_current_genbank"
        else:
            reason = "current_genbank_singleton"
        return selected, reason
    return None, "no_current_member"


def _optional_float(value: Any) -> float | None:
    if value in (None, ""):
        return None
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise MetadataQCError(f"Expected a numeric value, received {value!r}.") from exc
    if not math.isfinite(number):
        raise MetadataQCError(f"Numeric value must be finite, received {value!r}.")
    return number


def _optional_int(value: Any) -> int | None:
    number = _optional_float(value)
    if number is None:
        return None
    if not number.is_integer():
        raise MetadataQCError(f"Expected an integer value, received {value!r}.")
    return int(number)


def extract_selected_values(record: Mapping[str, Any]) -> dict[str, Any]:
    """Extract metadata fields from the newly selected current record."""
    organism = record.get("organism", {}) or {}
    assembly_info = record.get("assembly_info", {}) or {}
    assembly_stats = record.get("assembly_stats", {}) or {}
    return {
        "selected_accession": str(record["accession"]),
        "selected_source": record_source(record),
        "selected_status": record_status(record),
        "organism_name": organism.get("organism_name", ""),
        "tax_id": organism.get("tax_id", ""),
        "strain": (organism.get("infraspecific_names", {}) or {}).get("strain", ""),
        "assembly_level": assembly_info.get("assembly_level", ""),
        "genome_size_bp": assembly_stats.get("total_sequence_length", ""),
        "gc_percent": assembly_stats.get("gc_percent", ""),
        "number_of_contigs": assembly_stats.get("number_of_contigs", ""),
        "contig_n50": assembly_stats.get("contig_n50", ""),
    }


def modified_z_scores(values: Sequence[float | None]) -> tuple[list[float | None], float | None, float | None]:
    """Calculate modified z-scores without treating missing values as observations."""
    observed = [value for value in values if value is not None]
    if not observed:
        return [None] * len(values), None, None
    median = float(statistics.median(observed))
    mad = float(statistics.median(abs(value - median) for value in observed))
    if mad == 0:
        scores = [0.0 if value == median else None for value in values]
        return scores, median, mad
    scores = [
        None if value is None else 0.6745 * (value - median) / mad
        for value in values
    ]
    return scores, median, mad


def build_candidate_rows(
    summary_rows: Sequence[Mapping[str, str]],
    raw_records: Mapping[str, Mapping[str, Any]],
    current_statuses: set[str],
) -> list[dict[str, Any]]:
    """Resolve one available representative for each canonical assembly."""
    candidates: list[dict[str, Any]] = []
    for summary_row in summary_rows:
        original_accession = summary_row["representative_accession"]
        original = raw_records.get(original_accession)
        if original is None:
            raise MetadataQCError(
                f"Original representative {original_accession!r} is absent from raw metadata."
            )
        selected, reason = select_representative(
            summary_row, raw_records, current_statuses
        )
        base: dict[str, Any] = {
            "canonical_id": summary_row["canonical_id"],
            "gcf_accession": summary_row["gcf_accession"],
            "gca_accession": summary_row["gca_accession"],
            "pair_status": summary_row["pair_status"],
            "missing_partner_accession": summary_row["missing_partner_accession"],
            "original_representative_accession": original_accession,
            "original_representative_status": record_status(original),
            "selection_reason": reason,
        }
        if selected is None:
            base.update({field: "" for field in AUDIT_FIELDS if field not in base})
            base["metadata_decision"] = "exclude_unavailable"
            base["pangenome_contiguity_class"] = "unavailable"
            base["synteny_class"] = "unavailable"
            base["sequence_qc_status"] = "not_applicable"
            base["taxonomy_qc_status"] = "not_applicable"
            base["flags"] = "no_current_member"
        else:
            base.update(extract_selected_values(selected))
        candidates.append(base)
    return candidates


def classify_candidate(
    row: dict[str, Any],
    *,
    size_score: float | None,
    gc_score: float | None,
    threshold: float,
    contiguity: Mapping[str, Any],
) -> None:
    """Assign metadata, pangenome, and synteny classes to one available row."""
    if not row.get("selected_accession"):
        return
    flags: list[str] = []
    if row["selection_reason"] == "fallback_current_genbank":
        flags.append("suppressed_or_unavailable_refseq_fallback")
    if row["pair_status"] == "broken_pair":
        flags.append("broken_pair")
    required = {
        "organism_name": row.get("organism_name"),
        "tax_id": row.get("tax_id"),
        "genome_size_bp": row.get("genome_size_bp"),
        "number_of_contigs": row.get("number_of_contigs"),
        "contig_n50": row.get("contig_n50"),
    }
    flags.extend(
        f"missing_{field}" for field, value in required.items() if value in (None, "")
    )
    if row.get("gc_percent") in (None, ""):
        flags.append("missing_gc_percent")
    row["genome_size_modified_z"] = "" if size_score is None else round(size_score, 6)
    row["gc_percent_modified_z"] = "" if gc_score is None else round(gc_score, 6)
    if size_score is not None and abs(size_score) > threshold:
        flags.append("genome_size_outlier")
    if gc_score is not None and abs(gc_score) > threshold:
        flags.append("gc_percent_outlier")

    contigs = _optional_int(row.get("number_of_contigs"))
    n50 = _optional_int(row.get("contig_n50"))
    if contigs is None or n50 is None:
        row["pangenome_contiguity_class"] = "pending_metadata"
        row["synteny_class"] = "pending_metadata"
    else:
        pangenome_pass = (
            contigs <= int(contiguity["pangenome_max_contigs"])
            and n50 >= int(contiguity["pangenome_min_contig_n50_bp"])
        )
        synteny_pass = (
            contigs <= int(contiguity["synteny_max_contigs"])
            and n50 >= int(contiguity["synteny_min_contig_n50_bp"])
        )
        row["pangenome_contiguity_class"] = (
            "primary_contiguity_candidate" if pangenome_pass else "secondary_fragmented"
        )
        row["synteny_class"] = (
            "synteny_candidate" if synteny_pass else "presence_absence_only"
        )
        if not pangenome_pass:
            flags.append("pangenome_fragmented")
        if not synteny_pass:
            flags.append("synteny_fragmented")

    manual_review = any(
        flag.startswith("missing_") or flag.endswith("_outlier") for flag in flags
    )
    row["metadata_decision"] = "manual_review" if manual_review else "retain_candidate"
    row["sequence_qc_status"] = "pending"
    row["taxonomy_qc_status"] = "pending"
    row["flags"] = ";".join(sorted(set(flags)))


def build_report(
    rows: Sequence[Mapping[str, Any]],
    *,
    policy: Mapping[str, Any],
    raw_record_count: int,
    size_median: float | None,
    size_mad: float | None,
    gc_median: float | None,
    gc_mad: float | None,
) -> dict[str, Any]:
    """Summarize metadata-only decisions without claiming final sequence QC."""
    def counts(field: str) -> dict[str, int]:
        return dict(sorted(Counter(str(row[field]) for row in rows).items()))

    flag_counts = Counter(
        flag
        for row in rows
        for flag in str(row["flags"]).split(";")
        if flag
    )
    return {
        "schema_version": 1,
        "policy_id": policy["policy_id"],
        "audit_scope": "metadata_only",
        "final_qc_completed": False,
        "input_counts": {
            "canonical_assemblies": len(rows),
            "raw_accession_records": raw_record_count,
        },
        "selection_reason_counts": counts("selection_reason"),
        "metadata_decision_counts": counts("metadata_decision"),
        "pangenome_contiguity_counts": counts("pangenome_contiguity_class"),
        "synteny_class_counts": counts("synteny_class"),
        "flag_counts": dict(sorted(flag_counts.items())),
        "robust_baselines": {
            "genome_size_bp": {"median": size_median, "mad": size_mad},
            "gc_percent": {"median": gc_median, "mad": gc_mad},
            "modified_z_threshold": float(
                policy["robust_outlier_detection"]["threshold"]
            ),
        },
        "pending_steps": ["CheckM2", "GUNC", "GTDB-Tk_and_ANI_review"],
    }


def audit_metadata(
    summary_rows: Sequence[Mapping[str, str]],
    raw_records: Mapping[str, Mapping[str, Any]],
    policy: Mapping[str, Any],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Resolve representatives, flag anomalies, and classify metadata candidates."""
    validate_member_coverage(summary_rows, raw_records)
    current_statuses = {str(value).lower() for value in policy["current_status_values"]}
    base_rows = build_candidate_rows(summary_rows, raw_records, current_statuses)

    size_values = [
        _optional_float(row.get("genome_size_bp")) if row.get("selected_accession") else None
        for row in base_rows
    ]
    gc_values = [
        _optional_float(row.get("gc_percent")) if row.get("selected_accession") else None
        for row in base_rows
    ]
    size_scores, size_median, size_mad = modified_z_scores(size_values)
    gc_scores, gc_median, gc_mad = modified_z_scores(gc_values)
    threshold = float(policy["robust_outlier_detection"]["threshold"])
    contiguity = policy["metadata_contiguity"]

    for row, size_score, gc_score in zip(base_rows, size_scores, gc_scores, strict=True):
        classify_candidate(
            row,
            size_score=size_score,
            gc_score=gc_score,
            threshold=threshold,
            contiguity=contiguity,
        )
    summary = build_report(
        base_rows,
        policy=policy,
        raw_record_count=len(raw_records),
        size_median=size_median,
        size_mad=size_mad,
        gc_median=gc_median,
        gc_mad=gc_mad,
    )
    return base_rows, summary
