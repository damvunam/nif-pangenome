"""Regression tests for retrieval-manifest reconstruction."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from scripts import write_retrieval_manifest as manifest


class WriteRetrievalManifestTests(unittest.TestCase):
    """Verify computed counts, provenance preservation, and overwrite safety."""

    def test_build_manifest_preserves_retrieval_provenance(self) -> None:
        records = [
            {
                "accession": "GCA_000000001.1",
                "source_database": "SOURCE_DATABASE_GENBANK",
                "paired_accession": "GCF_000000001.1",
            },
            {
                "accession": "GCF_000000001.1",
                "source_database": "SOURCE_DATABASE_REFSEQ",
                "paired_accession": "GCA_000000001.1",
            },
            {
                "accession": "GCA_000000002.1",
                "source_database": "SOURCE_DATABASE_GENBANK",
            },
        ]
        provenance = {"tool": "NCBI Datasets CLI", "retrieval_timestamp": "fixed"}
        with tempfile.TemporaryDirectory() as directory:
            raw_path = Path(directory) / "raw.jsonl"
            raw_path.write_text("snapshot", encoding="utf-8")
            result = manifest.build_manifest(raw_path, records, provenance)

        self.assertEqual(provenance, result["retrieval"])
        self.assertEqual(3, result["empirical_counts"]["total_primary_records"])
        self.assertEqual(
            2, result["empirical_counts"]["records_with_paired_accession"]
        )
        self.assertEqual(
            1, result["empirical_counts"]["records_without_paired_accession"]
        )

    def test_duplicate_accession_is_rejected(self) -> None:
        item = {"accession": "GCA_000000001.1"}
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "raw.jsonl"
            path.write_text(
                json.dumps(item) + "\n" + json.dumps(item) + "\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(manifest.ManifestError, "Duplicate"):
                manifest.load_records(path)

    def test_writer_refuses_silent_overwrite(self) -> None:
        payload = {"retrieval": {"tool": "test"}}
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "manifest.json"
            path.write_text("sentinel", encoding="utf-8")
            with self.assertRaises(FileExistsError):
                manifest.write_manifest(payload, path, overwrite=False)
            self.assertEqual("sentinel", path.read_text(encoding="utf-8"))

            manifest.write_manifest(payload, path, overwrite=True)
            self.assertEqual(payload, json.loads(path.read_text(encoding="utf-8")))


if __name__ == "__main__":
    unittest.main()