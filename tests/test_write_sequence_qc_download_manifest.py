"""Tests for the sequence-QC pilot download manifest."""

from __future__ import annotations

import csv
import hashlib
import importlib.util
import tempfile
import unittest
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = REPOSITORY_ROOT / "scripts/write_sequence_qc_download_manifest.py"
COMMITTED_PILOT_PATH = (
    REPOSITORY_ROOT
    / "data/metadata/bradyrhizobium_sequence_qc_pilot_10.tsv"
)
COMMITTED_MANIFEST_PATH = (
    REPOSITORY_ROOT / "data/metadata/bradyrhizobium_download_manifest.tsv"
)

SPEC = importlib.util.spec_from_file_location(
    "write_sequence_qc_download_manifest", SCRIPT_PATH
)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"Cannot import {SCRIPT_PATH}.")
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def digest(path: Path, algorithm: str = "sha256") -> str:
    value = hashlib.new(algorithm)
    value.update(path.read_bytes())
    return value.hexdigest()


class DownloadManifestTests(unittest.TestCase):
    """Verify manifest validation, provenance, and committed reconciliation."""

    def make_fixture(self, root: Path) -> tuple[Path, Path]:
        pilot = root / "pilot.tsv"
        raw = root / "raw"
        data_root = raw / "unpacked/ncbi_dataset/data"
        data_root.mkdir(parents=True)

        pilot_rows = []
        accessions = []
        for index in range(1, 11):
            accession = f"GCF_{index:09d}.1"
            accessions.append(accession)
            pilot_rows.append(
                {
                    "pilot_order": index,
                    "pilot_role": f"role_{index}",
                    "canonical_id": f"CANON_{index:05d}",
                    "selected_accession": accession,
                }
            )
        with pilot.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=[
                    "pilot_order",
                    "pilot_role",
                    "canonical_id",
                    "selected_accession",
                ],
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(pilot_rows)

        (raw / "accessions.txt").write_text(
            "\n".join(accessions) + "\n", encoding="utf-8"
        )
        (raw / "datasets_version.txt").write_text(
            "datasets version: 18.33.1\n", encoding="utf-8"
        )
        (raw / "retrieval_timestamp_local.txt").write_text(
            "2026-08-15T12:00:00+07:00\n", encoding="utf-8"
        )
        (raw / "source_git_commit.txt").write_text(
            "a" * 40 + "\n", encoding="utf-8"
        )
        (raw / "source_pilot.sha256").write_text(
            f"{digest(pilot)}  pilot.tsv\n", encoding="utf-8"
        )

        package = raw / "ncbi_dataset.zip"
        package.write_bytes(b"synthetic package")
        (raw / "ncbi_dataset.zip.sha256").write_text(
            f"{digest(package)}  ncbi_dataset.zip\n", encoding="utf-8"
        )

        catalog = data_root / "dataset_catalog.json"
        catalog.write_text("{}\n", encoding="utf-8")
        md5_entries = [(digest(catalog, "md5"), "ncbi_dataset/data/dataset_catalog.json")]
        fasta_sha256_lines = []
        for accession in accessions:
            accession_root = data_root / accession
            accession_root.mkdir()
            fasta = accession_root / f"{accession}_genomic.fna"
            fasta.write_text(f">{accession}\nACGT\n", encoding="utf-8")
            package_path = f"ncbi_dataset/data/{accession}/{fasta.name}"
            md5_entries.append((digest(fasta, "md5"), package_path))
            fasta_sha256_lines.append(f"{digest(fasta)}  {package_path}")

        (raw / "unpacked/md5sum.txt").write_text(
            "\n".join(f"{checksum}  {path}" for checksum, path in md5_entries)
            + "\n",
            encoding="utf-8",
        )
        (raw / "genome_fasta_sha256.txt").write_text(
            "\n".join(fasta_sha256_lines) + "\n", encoding="utf-8"
        )
        return pilot, raw

    def test_build_manifest_validates_and_preserves_provenance(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            pilot, raw = self.make_fixture(Path(directory))
            rows = MODULE.build_download_manifest(
                MODULE.load_pilot(pilot), pilot, raw
            )
            self.assertEqual(10, len(rows))
            self.assertEqual("18.33.1", rows[0]["datasets_version"])
            self.assertEqual("pass", rows[0]["ncbi_md5_validation"])
            self.assertEqual(digest(pilot), rows[0]["source_pilot_sha256"])
            self.assertTrue(str(rows[0]["fasta_package_path"]).startswith("ncbi_dataset/"))

    def test_tampered_fasta_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            pilot, raw = self.make_fixture(Path(directory))
            fasta = next((raw / "unpacked/ncbi_dataset/data").glob("*/*_genomic.fna"))
            fasta.write_text(">tampered\nAAAA\n", encoding="utf-8")
            with self.assertRaisesRegex(MODULE.DownloadManifestError, "MD5 mismatch"):
                MODULE.build_download_manifest(MODULE.load_pilot(pilot), pilot, raw)

    def test_writer_refuses_silent_overwrite(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            pilot, raw = self.make_fixture(root)
            rows = MODULE.build_download_manifest(
                MODULE.load_pilot(pilot), pilot, raw
            )
            output = root / "manifest.tsv"
            MODULE.write_manifest_tsv(rows, output, overwrite=False)
            original = output.read_bytes()
            with self.assertRaises(FileExistsError):
                MODULE.write_manifest_tsv(rows, output, overwrite=False)
            self.assertEqual(original, output.read_bytes())

    def test_committed_manifest_reconciles_with_pilot_when_present(self) -> None:
        if not COMMITTED_MANIFEST_PATH.exists():
            self.skipTest("Download manifest has not been generated yet.")
        pilot_rows = MODULE.load_pilot(COMMITTED_PILOT_PATH)
        with COMMITTED_MANIFEST_PATH.open(newline="", encoding="utf-8") as handle:
            manifest_rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(10, len(manifest_rows))
        self.assertEqual(
            [row["selected_accession"] for row in pilot_rows],
            [row["selected_accession"] for row in manifest_rows],
        )
        self.assertEqual(
            [row["pilot_role"] for row in pilot_rows],
            [row["pilot_role"] for row in manifest_rows],
        )
        self.assertTrue(
            all(
                row["source_pilot_sha256"] == MODULE.file_digest(COMMITTED_PILOT_PATH)
                for row in manifest_rows
            )
        )


if __name__ == "__main__":
    unittest.main()
