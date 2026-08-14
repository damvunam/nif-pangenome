"""Unit tests for metadata-only biological QC."""

from __future__ import annotations

import csv
import json
import tempfile
import unittest
from pathlib import Path

from scripts import qc_metadata as qc


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
POLICY = qc.load_policy(REPOSITORY_ROOT / "config/metadata_qc_policy.json")


def raw_record(
    accession: str,
    *,
    status: str = "current",
    paired_accession: str | None = None,
    genome_size: int = 8_000_000,
    gc_percent: float | None = 63.5,
    contigs: int = 50,
    contig_n50: int = 100_000,
) -> dict[str, object]:
    """Create a compact NCBI-like assembly record."""
    source = (
        "SOURCE_DATABASE_REFSEQ"
        if accession.startswith("GCF_")
        else "SOURCE_DATABASE_GENBANK"
    )
    return {
        "accession": accession,
        "source_database": source,
        "paired_accession": paired_accession,
        "organism": {
            "organism_name": "Bradyrhizobium test",
            "tax_id": 374,
            "infraspecific_names": {"strain": accession},
        },
        "assembly_info": {
            "assembly_level": "Contig",
            "assembly_status": status,
        },
        "assembly_stats": {
            "total_sequence_length": genome_size,
            "gc_percent": gc_percent,
            "number_of_contigs": contigs,
            "contig_n50": contig_n50,
        },
    }


def summary_row(
    canonical_id: str,
    *,
    representative: str,
    gcf: str = "",
    gca: str = "",
    pair_status: str = "unpaired",
    missing_partner: str = "",
) -> dict[str, str]:
    """Create the summary fields required by metadata QC."""
    return {
        "canonical_id": canonical_id,
        "representative_accession": representative,
        "gcf_accession": gcf,
        "gca_accession": gca,
        "pair_status": pair_status,
        "missing_partner_accession": missing_partner,
    }


class RepresentativeSelectionTests(unittest.TestCase):
    """Verify current-member resolution and traceability."""

    def test_current_refseq_is_preferred(self) -> None:
        gcf = raw_record(
            "GCF_000000001.1", paired_accession="GCA_000000001.1"
        )
        gca = raw_record(
            "GCA_000000001.1", paired_accession="GCF_000000001.1"
        )
        row = summary_row(
            "CANON_00001",
            representative="GCF_000000001.1",
            gcf="GCF_000000001.1",
            gca="GCA_000000001.1",
            pair_status="paired",
        )
        selected, reason = qc.select_representative(
            row,
            {gcf["accession"]: gcf, gca["accession"]: gca},
            {"current", "latest"},
        )
        self.assertEqual("GCF_000000001.1", selected["accession"])
        self.assertEqual("current_refseq_preferred", reason)

    def test_suppressed_refseq_falls_back_to_current_genbank(self) -> None:
        gcf = raw_record(
            "GCF_000000001.1",
            status="suppressed",
            paired_accession="GCA_000000001.1",
        )
        gca = raw_record(
            "GCA_000000001.1", paired_accession="GCF_000000001.1"
        )
        row = summary_row(
            "CANON_00001",
            representative="GCF_000000001.1",
            gcf="GCF_000000001.1",
            gca="GCA_000000001.1",
            pair_status="paired",
        )
        selected, reason = qc.select_representative(
            row,
            {gcf["accession"]: gcf, gca["accession"]: gca},
            {"current", "latest"},
        )
        self.assertEqual("GCA_000000001.1", selected["accession"])
        self.assertEqual("fallback_current_genbank", reason)

    def test_no_current_member_is_unavailable(self) -> None:
        gca = raw_record("GCA_000000001.1", status="suppressed")
        row = summary_row(
            "CANON_00001",
            representative="GCA_000000001.1",
            gca="GCA_000000001.1",
        )
        selected, reason = qc.select_representative(
            row, {gca["accession"]: gca}, {"current", "latest"}
        )
        self.assertIsNone(selected)
        self.assertEqual("no_current_member", reason)


class MetadataAuditTests(unittest.TestCase):
    """Verify cohort classification and deterministic validation."""

    def test_audit_resolves_fallback_and_broken_pair(self) -> None:
        paired_gcf = raw_record(
            "GCF_000000001.1", paired_accession="GCA_000000001.1"
        )
        paired_gca = raw_record(
            "GCA_000000001.1", paired_accession="GCF_000000001.1"
        )
        suppressed_gcf = raw_record(
            "GCF_000000002.1",
            status="suppressed",
            paired_accession="GCA_000000002.1",
        )
        fallback_gca = raw_record(
            "GCA_000000002.1",
            paired_accession="GCF_000000002.1",
            contigs=600,
            contig_n50=10_000,
        )
        broken_gca = raw_record(
            "GCA_000000003.1", paired_accession="GCF_000000003.1"
        )
        records = {
            str(record["accession"]): record
            for record in (
                paired_gcf,
                paired_gca,
                suppressed_gcf,
                fallback_gca,
                broken_gca,
            )
        }
        rows = [
            summary_row(
                "CANON_00001",
                representative="GCF_000000001.1",
                gcf="GCF_000000001.1",
                gca="GCA_000000001.1",
                pair_status="paired",
            ),
            summary_row(
                "CANON_00002",
                representative="GCF_000000002.1",
                gcf="GCF_000000002.1",
                gca="GCA_000000002.1",
                pair_status="paired",
            ),
            summary_row(
                "CANON_00003",
                representative="GCA_000000003.1",
                gca="GCA_000000003.1",
                pair_status="broken_pair",
                missing_partner="GCF_000000003.1",
            ),
        ]

        audit, report = qc.audit_metadata(rows, records, POLICY)

        self.assertEqual("GCA_000000002.1", audit[1]["selected_accession"])
        self.assertIn("suppressed_or_unavailable_refseq_fallback", audit[1]["flags"])
        self.assertEqual("secondary_fragmented", audit[1]["pangenome_contiguity_class"])
        self.assertEqual("presence_absence_only", audit[1]["synteny_class"])
        self.assertIn("broken_pair", audit[2]["flags"])
        self.assertEqual(3, report["input_counts"]["canonical_assemblies"])
        self.assertFalse(report["final_qc_completed"])

    def test_missing_required_metadata_requires_manual_review(self) -> None:
        gca = raw_record("GCA_000000001.1")
        gca["assembly_stats"]["contig_n50"] = None
        row = summary_row(
            "CANON_00001",
            representative="GCA_000000001.1",
            gca="GCA_000000001.1",
        )
        audit, _ = qc.audit_metadata([row], {gca["accession"]: gca}, POLICY)
        self.assertEqual("manual_review", audit[0]["metadata_decision"])
        self.assertIn("missing_contig_n50", audit[0]["flags"])

    def test_unavailable_member_is_excluded_without_final_qc_claim(self) -> None:
        gca = raw_record("GCA_000000001.1", status="suppressed")
        row = summary_row(
            "CANON_00001",
            representative="GCA_000000001.1",
            gca="GCA_000000001.1",
        )
        audit, report = qc.audit_metadata([row], {gca["accession"]: gca}, POLICY)
        self.assertEqual("exclude_unavailable", audit[0]["metadata_decision"])
        self.assertEqual("not_applicable", audit[0]["sequence_qc_status"])
        self.assertFalse(report["final_qc_completed"])

    def test_raw_and_summary_accessions_must_match_exactly(self) -> None:
        represented = raw_record("GCA_000000001.1")
        extra = raw_record("GCA_000000002.1")
        row = summary_row(
            "CANON_00001",
            representative="GCA_000000001.1",
            gca="GCA_000000001.1",
        )
        with self.assertRaisesRegex(qc.MetadataQCError, "absent from the canonical"):
            qc.audit_metadata(
                [row],
                {represented["accession"]: represented, extra["accession"]: extra},
                POLICY,
            )

    def test_modified_z_scores_are_deterministic(self) -> None:
        first = qc.modified_z_scores([1.0, 2.0, 3.0, 4.0, 100.0, None])
        second = qc.modified_z_scores([1.0, 2.0, 3.0, 4.0, 100.0, None])
        self.assertEqual(first, second)
        self.assertGreater(abs(first[0][4]), 3.5)
        self.assertIsNone(first[0][5])

    def test_atomic_writers_refuse_silent_overwrite(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            tsv = Path(directory) / "audit.tsv"
            report = Path(directory) / "report.json"
            tsv.write_text("sentinel", encoding="utf-8")
            report.write_text("sentinel", encoding="utf-8")
            with self.assertRaises(FileExistsError):
                qc.write_tsv_atomic([], tsv, overwrite=False)
            with self.assertRaises(FileExistsError):
                qc.write_json_atomic({}, report, overwrite=False)

    def test_audit_tsv_uses_declared_schema(self) -> None:
        gca = raw_record("GCA_000000001.1")
        row = summary_row(
            "CANON_00001",
            representative="GCA_000000001.1",
            gca="GCA_000000001.1",
        )
        audit, _ = qc.audit_metadata([row], {gca["accession"]: gca}, POLICY)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "audit.tsv"
            qc.write_tsv_atomic(audit, path, overwrite=False)
            with path.open(newline="", encoding="utf-8") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                written = list(reader)
            self.assertEqual(qc.AUDIT_FIELDS, reader.fieldnames)
            self.assertEqual(1, len(written))


if __name__ == "__main__":
    unittest.main()