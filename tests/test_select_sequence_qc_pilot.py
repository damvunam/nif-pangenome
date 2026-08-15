"""Tests for deterministic sequence-QC pilot selection."""

from __future__ import annotations

import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = REPOSITORY_ROOT / "scripts/select_sequence_qc_pilot.py"
INPUT_PATH = REPOSITORY_ROOT / "data/metadata/bradyrhizobium_metadata_qc.tsv"
POLICY_PATH = REPOSITORY_ROOT / "config/metadata_qc_policy.json"
COMMITTED_PILOT_PATH = (
    REPOSITORY_ROOT
    / "data/metadata/bradyrhizobium_sequence_qc_pilot_10.tsv"
)

SPEC = importlib.util.spec_from_file_location("select_sequence_qc_pilot", SCRIPT_PATH)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"Cannot import {SCRIPT_PATH}.")
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class SequenceQCPilotTests(unittest.TestCase):
    """Verify deterministic strata, committed output, and safe writing."""

    @classmethod
    def setUpClass(cls) -> None:
        cls.rows = MODULE.load_metadata_qc(INPUT_PATH)
        cls.policy = MODULE.load_contiguity_policy(POLICY_PATH)
        cls.selected = MODULE.select_pilot_rows(cls.rows, cls.policy)
        cls.output_rows = MODULE.build_output_rows(
            cls.selected,
            metadata_qc_sha256=MODULE.sha256_file(INPUT_PATH),
            policy_sha256=MODULE.sha256_file(POLICY_PATH),
        )

    def test_pilot_has_ten_unique_accessions_and_locked_roles(self) -> None:
        expected_roles = [
            "normal_refseq_near_median_1",
            "normal_refseq_near_median_2",
            "normal_genbank_singleton_near_median",
            "fragmented_extreme",
            "genome_size_outlier_isolated",
            "gc_percent_outlier_isolated",
            "missing_gc_percent_isolated",
            "refseq_to_genbank_fallback_isolated",
            "broken_pair",
            "synteny_threshold_edge",
        ]
        self.assertEqual(expected_roles, [row["pilot_role"] for row in self.output_rows])
        accessions = [row["selected_accession"] for row in self.output_rows]
        self.assertEqual(10, len(accessions))
        self.assertEqual(10, len(set(accessions)))

    def test_selection_is_independent_of_input_order(self) -> None:
        reversed_selected = MODULE.select_pilot_rows(
            list(reversed(self.rows)), self.policy
        )
        self.assertEqual(
            [row[2]["selected_accession"] for row in self.selected],
            [row[2]["selected_accession"] for row in reversed_selected],
        )

    def test_committed_pilot_matches_current_inputs(self) -> None:
        with COMMITTED_PILOT_PATH.open(newline="", encoding="utf-8") as handle:
            committed = list(csv.DictReader(handle, delimiter="\t"))
        expected = [
            {field: str(row[field]) for field in MODULE.PILOT_FIELDS}
            for row in self.output_rows
        ]
        self.assertEqual(expected, committed)

    def test_missing_stratum_is_rejected(self) -> None:
        without_broken_pairs = [
            row for row in self.rows if "broken_pair" not in row["flags"].split(";")
        ]
        with self.assertRaisesRegex(MODULE.PilotSelectionError, "broken_pair"):
            MODULE.select_pilot_rows(without_broken_pairs, self.policy)

    def test_writer_refuses_silent_overwrite(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "pilot.tsv"
            MODULE.write_pilot_tsv(self.output_rows, output, overwrite=False)
            original = output.read_bytes()
            with self.assertRaises(FileExistsError):
                MODULE.write_pilot_tsv(self.output_rows, output, overwrite=False)
            self.assertEqual(original, output.read_bytes())


if __name__ == "__main__":
    unittest.main()
