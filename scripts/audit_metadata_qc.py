"""Run metadata-only biological QC for canonical Bradyrhizobium assemblies."""

from __future__ import annotations

import argparse
from collections.abc import Sequence
from pathlib import Path

from qc_metadata import (
    audit_metadata,
    load_policy,
    load_raw_records,
    load_summary,
    sha256_file,
    write_json_atomic,
    write_tsv_atomic,
)


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_RAW_PATH = REPOSITORY_ROOT / "data/raw/bradyrhizobium_genome_summary.jsonl"
DEFAULT_SUMMARY_PATH = (
    REPOSITORY_ROOT / "data/metadata/bradyrhizobium_assembly_summary.tsv"
)
DEFAULT_POLICY_PATH = REPOSITORY_ROOT / "config/metadata_qc_policy.json"
DEFAULT_AUDIT_OUTPUT = (
    REPOSITORY_ROOT / "data/metadata/bradyrhizobium_metadata_qc.tsv"
)
DEFAULT_REPORT_OUTPUT = REPOSITORY_ROOT / "reports/task2_metadata_qc_summary.json"


def preflight_outputs(paths: Sequence[Path], *, overwrite: bool) -> None:
    """Reject conflicting outputs before either artifact is written."""
    if overwrite:
        return
    existing = [str(path) for path in paths if path.exists()]
    if existing:
        raise FileExistsError(
            f"Output artifacts already exist: {existing!r}. Use --overwrite explicitly."
        )


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-path", type=Path, default=DEFAULT_RAW_PATH)
    parser.add_argument("--summary-path", type=Path, default=DEFAULT_SUMMARY_PATH)
    parser.add_argument("--policy-path", type=Path, default=DEFAULT_POLICY_PATH)
    parser.add_argument("--audit-output", type=Path, default=DEFAULT_AUDIT_OUTPUT)
    parser.add_argument("--report-output", type=Path, default=DEFAULT_REPORT_OUTPUT)
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Explicitly replace both output artifacts atomically.",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    """Run the audit, save deterministic artifacts, and print decision counts."""
    args = parse_args(argv)
    policy = load_policy(args.policy_path)
    summary_rows = load_summary(args.summary_path)
    raw_records = load_raw_records(args.raw_path)
    audit_rows, report = audit_metadata(summary_rows, raw_records, policy)
    report["input_fingerprints"] = {
        "raw_sha256": sha256_file(args.raw_path),
        "assembly_summary_sha256": sha256_file(args.summary_path),
        "policy_sha256": sha256_file(args.policy_path),
    }
    preflight_outputs(
        [args.audit_output, args.report_output], overwrite=args.overwrite
    )
    write_tsv_atomic(audit_rows, args.audit_output, overwrite=args.overwrite)
    write_json_atomic(report, args.report_output, overwrite=args.overwrite)
    print(f"Metadata QC rows: {len(audit_rows)}")
    print("Selection reasons:", report["selection_reason_counts"])
    print("Metadata decisions:", report["metadata_decision_counts"])
    print("Pangenome classes:", report["pangenome_contiguity_counts"])
    print("Synteny classes:", report["synteny_class_counts"])
    print("Final sequence QC completed: NO")
    print(f"Audit TSV: {args.audit_output}")
    print(f"Summary JSON: {args.report_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
