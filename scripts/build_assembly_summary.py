"""
Builds data/metadata/bradyrhizobium_assembly_summary.tsv — one row per
CANONICAL PHYSICAL ASSEMBLY (GCA/GCF pairs collapsed to a single row).

Representative selection rule (per confirmed decision):
  - If a GCA/GCF pair exists, use the GCF (RefSeq) record as the primary
    source for descriptive/quality fields, and record the paired GCA
    accession alongside for traceability.
  - If no pair exists (singleton group), use whichever single record
    is present as-is.

Input: data/raw/bradyrhizobium_genome_summary.jsonl (gitignored, local only)
Output: data/metadata/bradyrhizobium_assembly_summary.tsv (git-tracked)

This script performs NO filtering, NO QC thresholding, and NO exclusion
of any record. All 2211 canonical groups present in the raw data are
written as rows, including suppressed-status and broken-pair records,
which are flagged via explicit columns rather than dropped silently.
"""
import json
import csv

RAW_PATH = "data/raw/bradyrhizobium_genome_summary.jsonl"
OUT_PATH = "data/metadata/bradyrhizobium_assembly_summary.tsv"

FIELDS = [
    "canonical_id",
    "representative_accession",
    "representative_source",       # REFSEQ or GENBANK
    "gcf_accession",
    "gca_accession",
    "pair_status",                 # paired | unpaired | broken_pair
    "organism_name",
    "tax_id",
    "strain",
    "assembly_level",
    "assembly_status",
    "refseq_category",
    "genome_size_bp",
    "gc_percent",
    "number_of_contigs",
    "number_of_scaffolds",
    "contig_n50",
    "scaffold_n50",
    "bioproject_accession",
    "biosample_accession",
    "submitter",
    "release_date",
]

def extract_row(rec, canonical_id, gcf_acc, gca_acc, pair_status):
    organism = rec.get("organism", {}) or {}
    assembly_info = rec.get("assembly_info", {}) or {}
    assembly_stats = rec.get("assembly_stats", {}) or {}
    biosample = assembly_info.get("biosample", {}) or {}

    return {
        "canonical_id": canonical_id,
        "representative_accession": rec.get("accession"),
        "representative_source": rec.get("source_database", "").replace("SOURCE_DATABASE_", ""),
        "gcf_accession": gcf_acc or "",
        "gca_accession": gca_acc or "",
        "pair_status": pair_status,
        "organism_name": organism.get("organism_name", ""),
        "tax_id": organism.get("tax_id", ""),
        "strain": (organism.get("infraspecific_names", {}) or {}).get("strain", ""),
        "assembly_level": assembly_info.get("assembly_level", ""),
        "assembly_status": assembly_info.get("assembly_status", ""),
        "refseq_category": assembly_info.get("refseq_category", ""),
        "genome_size_bp": assembly_stats.get("total_sequence_length", ""),
        "gc_percent": assembly_stats.get("gc_percent", ""),
        "number_of_contigs": assembly_stats.get("number_of_contigs", ""),
        "number_of_scaffolds": assembly_stats.get("number_of_scaffolds", ""),
        "contig_n50": assembly_stats.get("contig_n50", ""),
        "scaffold_n50": assembly_stats.get("scaffold_n50", ""),
        "bioproject_accession": assembly_info.get("bioproject_accession", ""),
        "biosample_accession": biosample.get("accession", ""),
        "submitter": assembly_info.get("submitter", ""),
        "release_date": assembly_info.get("release_date", ""),
    }

def main():
    records = {}
    with open(RAW_PATH) as f:
        for line in f:
            r = json.loads(line)
            records[r["accession"]] = r

    seen = set()
    rows = []
    canonical_counter = 0

    for acc, rec in records.items():
        if acc in seen:
            continue

        partner_acc = rec.get("paired_accession")
        partner_rec = records.get(partner_acc) if partner_acc else None

        if partner_acc and partner_rec is not None:
            # Genuine pair present in dataset
            group = {acc, partner_acc}
            pair_status = "paired"
            # Determine which is GCF vs GCA
            if rec.get("source_database") == "SOURCE_DATABASE_REFSEQ":
                refseq_rec, genbank_acc = rec, partner_acc
            else:
                refseq_rec, genbank_acc = partner_rec, acc
            gcf_acc = refseq_rec.get("accession")
            gca_acc = genbank_acc
            representative = refseq_rec

        elif partner_acc and partner_rec is None:
            # Broken pair: claims a partner not present in this dataset
            group = {acc}
            pair_status = "broken_pair"
            representative = rec
            if rec.get("source_database") == "SOURCE_DATABASE_REFSEQ":
                gcf_acc, gca_acc = rec.get("accession"), ""
            else:
                gcf_acc, gca_acc = "", rec.get("accession")

        else:
            # No paired_accession at all
            group = {acc}
            pair_status = "unpaired"
            representative = rec
            if rec.get("source_database") == "SOURCE_DATABASE_REFSEQ":
                gcf_acc, gca_acc = rec.get("accession"), ""
            else:
                gcf_acc, gca_acc = "", rec.get("accession")

        canonical_counter += 1
        canonical_id = f"CANON_{canonical_counter:05d}"

        rows.append(extract_row(representative, canonical_id, gcf_acc, gca_acc, pair_status))
        seen |= group

    with open(OUT_PATH, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} canonical assembly rows to {OUT_PATH}")
    print(f"Total input records consumed: {len(seen)} (should equal {len(records)})")

    # Quick pair_status breakdown for sanity check
    from collections import Counter
    print("pair_status breakdown:", Counter(r["pair_status"] for r in rows))
    print("representative_source breakdown:", Counter(r["representative_source"] for r in rows))

if __name__ == "__main__":
    main()
