"""Integrity tests for the metadata artifacts committed with Task 1."""

from __future__ import annotations

import csv
import json
import unittest
from collections import Counter
from pathlib import Path

from scripts import build_assembly_summary as summary


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
SUMMARY_PATH = REPOSITORY_ROOT / "data/metadata/bradyrhizobium_assembly_summary.tsv"
MANIFEST_PATH = REPOSITORY_ROOT / "data/metadata/retrieval_manifest.json"


class CommittedMetadataTests(unittest.TestCase):
    """Ensure checked-in metadata remains internally consistent."""

    @classmethod
    def setUpClass(cls) -> None:
        with SUMMARY_PATH.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            cls.rows = list(reader)
            cls.fields = reader.fieldnames
        cls.manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))

    def test_schema_and_canonical_counts(self) -> None:
        self.assertEqual(summary.FIELDS, self.fields)
        self.assertEqual(2_211, len(self.rows))
        self.assertEqual(2_211, len({row["canonical_id"] for row in self.rows}))
        self.assertEqual(
            [f"CANON_{index:05d}" for index in range(1, 2_212)],
            [row["canonical_id"] for row in self.rows],
        )

    def test_rows_use_stable_accession_order(self) -> None:
        keys = [row["gca_accession"] or row["gcf_accession"] for row in self.rows]
        self.assertEqual(sorted(keys), keys)

    def test_pair_breakdown_and_accession_uniqueness(self) -> None:
        self.assertEqual(
            Counter({"paired": 1_876, "unpaired": 331, "broken_pair": 4}),
            Counter(row["pair_status"] for row in self.rows),
        )
        gcf = [row["gcf_accession"] for row in self.rows if row["gcf_accession"]]
        gca = [row["gca_accession"] for row in self.rows if row["gca_accession"]]
        self.assertEqual(1_876, len(gcf))
        self.assertEqual(2_211, len(gca))
        self.assertEqual(len(gcf), len(set(gcf)))
        self.assertEqual(len(gca), len(set(gca)))

    def test_broken_pair_traceability(self) -> None:
        expected = {
            "GCF_059803835.1",
            "GCF_059803455.1",
            "GCF_059803515.1",
            "GCF_059803475.1",
        }
        broken = [row for row in self.rows if row["pair_status"] == "broken_pair"]
        self.assertEqual(expected, {row["missing_partner_accession"] for row in broken})
        self.assertTrue(all(not row["gcf_accession"] for row in broken))
        self.assertTrue(
            all(not row["missing_partner_accession"] for row in self.rows if row not in broken)
        )

    def test_manifest_reconciles_with_summary(self) -> None:
        counts = self.manifest["empirical_counts"]
        paired = sum(row["pair_status"] == "paired" for row in self.rows)
        broken = sum(row["pair_status"] == "broken_pair" for row in self.rows)
        unpaired = sum(row["pair_status"] == "unpaired" for row in self.rows)
        self.assertEqual(counts["total_primary_records"], paired * 2 + broken + unpaired)
        self.assertEqual(counts["records_with_paired_accession"], paired * 2 + broken)
        self.assertEqual(counts["records_without_paired_accession"], unpaired)


if __name__ == "__main__":
    unittest.main()