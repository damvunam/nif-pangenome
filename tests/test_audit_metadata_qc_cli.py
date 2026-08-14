"""Integration tests for the metadata-QC command-line interface."""

from __future__ import annotations

import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPOSITORY_ROOT / "scripts/audit_metadata_qc.py"
POLICY = REPOSITORY_ROOT / "config/metadata_qc_policy.json"


class MetadataQCCLITests(unittest.TestCase):
    """Verify successful generation and all-output overwrite protection."""

    def test_cli_generates_both_artifacts_and_refuses_partial_overwrite(self) -> None:
        record = {
            "accession": "GCA_000000001.1",
            "source_database": "SOURCE_DATABASE_GENBANK",
            "organism": {"organism_name": "Bradyrhizobium test", "tax_id": 374},
            "assembly_info": {
                "assembly_level": "Contig",
                "assembly_status": "current",
            },
            "assembly_stats": {
                "total_sequence_length": 8_000_000,
                "gc_percent": 63.5,
                "number_of_contigs": 20,
                "contig_n50": 200_000,
            },
        }
        summary_fields = [
            "canonical_id",
            "representative_accession",
            "gcf_accession",
            "gca_accession",
            "pair_status",
            "missing_partner_accession",
        ]
        summary_row = {
            "canonical_id": "CANON_00001",
            "representative_accession": "GCA_000000001.1",
            "gcf_accession": "",
            "gca_accession": "GCA_000000001.1",
            "pair_status": "unpaired",
            "missing_partner_accession": "",
        }

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            raw = root / "raw.jsonl"
            summary = root / "summary.tsv"
            audit = root / "audit.tsv"
            report = root / "report.json"
            raw.write_text(json.dumps(record) + "\n", encoding="utf-8")
            with summary.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(
                    handle, fieldnames=summary_fields, delimiter="\t", lineterminator="\n"
                )
                writer.writeheader()
                writer.writerow(summary_row)

            command = [
                sys.executable,
                str(SCRIPT),
                "--raw-path",
                str(raw),
                "--summary-path",
                str(summary),
                "--policy-path",
                str(POLICY),
                "--audit-output",
                str(audit),
                "--report-output",
                str(report),
            ]
            completed = subprocess.run(
                command,
                cwd=REPOSITORY_ROOT,
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertEqual(0, completed.returncode, completed.stderr)
            self.assertIn("Final sequence QC completed: NO", completed.stdout)
            self.assertTrue(audit.exists())
            payload = json.loads(report.read_text(encoding="utf-8"))
            self.assertFalse(payload["final_qc_completed"])

            audit.unlink()
            second = subprocess.run(
                command,
                cwd=REPOSITORY_ROOT,
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertNotEqual(0, second.returncode)
            self.assertFalse(audit.exists(), "Preflight must prevent a partial TSV write.")


if __name__ == "__main__":
    unittest.main()
