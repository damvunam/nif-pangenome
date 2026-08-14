"""Regression tests for deterministic canonical assembly summarization."""

from __future__ import annotations

import csv
import json
import tempfile
import unittest
from pathlib import Path

from scripts import build_assembly_summary as summary


def record(
    accession: str,
    source: str,
    paired_accession: str | None = None,
) -> dict[str, object]:
    """Create the smallest valid assembly record used by these tests."""
    return {
        "accession": accession,
        "source_database": source,
        "paired_accession": paired_accession,
        "organism": {"organism_name": "Bradyrhizobium test", "tax_id": 374},
        "assembly_info": {"assembly_level": "Contig", "assembly_status": "current"},
        "assembly_stats": {"total_sequence_length": 10, "number_of_contigs": 1},
    }


class BuildAssemblySummaryTests(unittest.TestCase):
    """Verify grouping, validation, deterministic order, and safe output."""

    def test_rows_are_stable_when_input_order_changes(self) -> None:
        gca = record(
            "GCA_000000002.1",
            summary.SOURCE_GENBANK,
            "GCF_000000002.1",
        )
        gcf = record(
            "GCF_000000002.1",
            summary.SOURCE_REFSEQ,
            "GCA_000000002.1",
        )
        singleton = record("GCA_000000001.1", summary.SOURCE_GENBANK)
        first = {item["accession"]: item for item in (gca, gcf, singleton)}
        second = {item["accession"]: item for item in (singleton, gcf, gca)}

        self.assertEqual(summary.build_rows(first), summary.build_rows(second))
        rows = summary.build_rows(first)
        self.assertEqual("GCA_000000001.1", rows[0]["gca_accession"])
        self.assertEqual("GCF_000000002.1", rows[1]["representative_accession"])

    def test_broken_pair_preserves_missing_partner_accession(self) -> None:
        gca = record(
            "GCA_059803835.1",
            summary.SOURCE_GENBANK,
            "GCF_059803835.1",
        )
        row = summary.build_rows({str(gca["accession"]): gca})[0]

        self.assertEqual("broken_pair", row["pair_status"])
        self.assertEqual("GCF_059803835.1", row["missing_partner_accession"])
        self.assertEqual("", row["gcf_accession"])

    def test_nonreciprocal_pair_is_rejected(self) -> None:
        gca = record(
            "GCA_000000001.1",
            summary.SOURCE_GENBANK,
            "GCF_000000001.1",
        )
        gcf = record(
            "GCF_000000001.1",
            summary.SOURCE_REFSEQ,
            "GCA_999999999.1",
        )
        records = {str(item["accession"]): item for item in (gca, gcf)}

        with self.assertRaisesRegex(summary.AssemblySummaryError, "Non-reciprocal"):
            summary.build_rows(records)

    def test_broken_pair_with_same_database_prefix_is_rejected(self) -> None:
        gca = record(
            "GCA_000000001.1",
            summary.SOURCE_GENBANK,
            "GCA_000000002.1",
        )

        with self.assertRaisesRegex(summary.AssemblySummaryError, "must start with"):
            summary.build_rows({str(gca["accession"]): gca})

    def test_duplicate_accession_is_rejected(self) -> None:
        duplicate = record("GCA_000000001.1", summary.SOURCE_GENBANK)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "raw.jsonl"
            path.write_text(
                json.dumps(duplicate) + "\n" + json.dumps(duplicate) + "\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(summary.AssemblySummaryError, "Duplicate"):
                summary.load_records(path)

    def test_writer_refuses_silent_overwrite(self) -> None:
        row = summary.build_rows(
            {
                "GCA_000000001.1": record(
                    "GCA_000000001.1", summary.SOURCE_GENBANK
                )
            }
        )[0]
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "summary.tsv"
            path.write_text("sentinel", encoding="utf-8")
            with self.assertRaises(FileExistsError):
                summary.write_rows([row], path, overwrite=False)
            self.assertEqual("sentinel", path.read_text(encoding="utf-8"))

            summary.write_rows([row], path, overwrite=True)
            with path.open(newline="", encoding="utf-8") as handle:
                written = list(csv.DictReader(handle, delimiter="\t"))
            expected = {key: str(value) for key, value in row.items()}
            self.assertEqual([expected], written)


if __name__ == "__main__":
    unittest.main()