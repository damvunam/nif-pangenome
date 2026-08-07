"""
Writes data/metadata/retrieval_manifest.json documenting the raw NCBI
genome summary retrieval for Bradyrhizobium (TaxID 374).

All values here are either:
  - the literal command executed,
  - directly measured from the retrieved file (size, checksum, line count),
  - or counts computed from the file's own content (record/pair statistics).
Nothing here is copied from prior/unverified reports.
"""
import json
import hashlib
import os

RAW_PATH = "data/raw/bradyrhizobium_genome_summary.jsonl"
MANIFEST_PATH = "data/metadata/retrieval_manifest.json"

def sha256sum(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

def main():
    size_bytes = os.path.getsize(RAW_PATH)
    checksum = sha256sum(RAW_PATH)

    with open(RAW_PATH) as f:
        lines = f.readlines()
    total_records = len(lines)

    records = [json.loads(l) for l in lines]
    accessions = [r.get("accession") for r in records]
    distinct_accessions = len(set(accessions))

    source_counts = {}
    for r in records:
        sd = r.get("source_database", "MISSING")
        source_counts[sd] = source_counts.get(sd, 0) + 1

    paired_count = sum(1 for r in records if r.get("paired_accession"))
    unpaired_count = total_records - paired_count

    manifest = {
        "retrieval": {
            "tool": "NCBI Datasets CLI",
            "command": (
                "datasets summary genome taxon 374 "
                "--as-json-lines --assembly-version current "
                "> data/raw/bradyrhizobium_genome_summary.jsonl"
            ),
            "tool_version_used": "18.33.1",
            "tool_version_available_at_time_of_run": "18.34.0",
            "taxon_scientific_name": "Bradyrhizobium",
            "taxon_id": 374,
            "assembly_version_scope": "current",
            "retrieved_by": "damvunam",
            "retrieval_timestamp_local": "2026-08-07T14:43 (Asia/Ho_Chi_Minh, unconfirmed timezone offset — file mtime based)",
            "note_on_timestamp": "Timestamp derived from file mtime reported by `ls -la` at retrieval time; not independently verified against system clock/timezone."
        },
        "raw_file": {
            "path": RAW_PATH,
            "size_bytes": size_bytes,
            "sha256": checksum,
            "git_tracked": False,
            "note": "Excluded from git via .gitignore; integrity tracked here via checksum."
        },
        "empirical_counts": {
            "total_primary_records": total_records,
            "distinct_accessions": distinct_accessions,
            "source_database_counts": source_counts,
            "records_with_paired_accession": paired_count,
            "records_without_paired_accession": unpaired_count
        }
    }

    with open(MANIFEST_PATH, "w") as f:
        json.dump(manifest, f, indent=2)

    print(f"Manifest written to {MANIFEST_PATH}")
    print(json.dumps(manifest, indent=2))

if __name__ == "__main__":
    main()
