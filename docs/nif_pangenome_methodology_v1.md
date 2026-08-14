# Nif-Pangenome Methodology v1

**Version:** 1.0.0
**Date:** 2026-08-15
**Status:** Approved design; sequence-based QC execution pending
**Primary biological scope:** Comparative pangenomics of `nif`, `nod`, and
`fix` systems across *Bradyrhizobium* genomes

## 1. Scientific objective

Nif-Pangenome investigates the distribution, conservation, structural
variation, and evolutionary context of nitrogen-fixation and symbiosis genes
across the *Bradyrhizobium* pangenome. Genome retrieval and quality control are
prerequisites rather than the final scientific objective.

The central questions are:

1. Which genes are core, accessory, or rare across *Bradyrhizobium* genomes?
2. How complete and structurally conserved are the `nif`, `nod`, and `fix`
   systems?
3. Which variations reflect biological evolution, and which may be caused by
   incomplete, fragmented, contaminated, or taxonomically inconsistent
   assemblies?
4. How do symbiosis-gene histories compare with whole-genome taxonomy and the
   core-genome phylogeny?

## 2. Method status vocabulary

Each important rule is identified as one of the following:

- **Adopted:** used according to a published method or official tool policy.
- **Adapted:** derived from a published standard but modified for this project.
- **Nif-Pangenome operational rule:** a project-specific safeguard that must
  not be presented as a universal biological threshold.

## 3. Current data scope and provenance

The metadata snapshot was retrieved from NCBI Taxonomy ID 374
(*Bradyrhizobium*) using accession-version identifiers.

The audited snapshot contains:

- 4,087 primary NCBI assembly records;
- 1,876 RefSeq records and 2,211 GenBank records;
- 2,211 canonical physical assemblies after reciprocal GCF/GCA pairs were
  collapsed;
- 1,876 reciprocal GCF/GCA pairs;
- 331 GenBank-only assemblies;
- 4 broken-pair records retained with explicit traceability;
- 58 suppressed RefSeq representatives that fall back to current GenBank
  partners.

Raw files are excluded from Git. Their paths, sizes, retrieval commands, tool
versions, timestamps, and SHA-256 checksums are retained in manifests.

For more than 1,000 genomes, sequence retrieval will follow the NCBI Datasets
dehydrated-download and rehydration workflow [1]. Exact accession versions must
be used. A suppressed or unavailable accession must never be silently replaced
during download.

## 4. Canonical assembly selection

This section contains **Nif-Pangenome operational rules**.

For each canonical physical assembly:

1. Prefer a current RefSeq (`GCF_`) member.
2. If RefSeq is suppressed but its GenBank partner is current, select the
   current `GCA_` record and retain the fallback reason.
3. Retain a current broken-pair GenBank assembly with an explicit traceability
   flag.
4. Mark an assembly unavailable only when no current GCF or GCA member remains.
5. Preserve the original representative accession, selected accession, member
   statuses, pair status, and selection reason.

The output must contain one unique selected accession for each available
canonical assembly.

## 5. Metadata-only screening

Metadata screening is not final genome QC.

### 5.1 Contiguity classes

The following thresholds are **Nif-Pangenome operational rules**:

- **Primary pangenome contiguity candidate:** no more than 500 contigs and
  contig N50 of at least 20 kb.
- **Synteny candidate:** no more than 200 contigs and contig N50 of at least
  50 kb.
- **Secondary fragmented candidate:** available but failing one or more of the
  relevant contiguity thresholds.

N50 is not treated as a standalone measure of assembly correctness because it
can improve despite misassembly [2]. Fragmentation may also inflate apparent
gene absence and distort core/accessory genome estimates [3].

### 5.2 Robust metadata outlier screen

Genome size and GC percentage are screened with the **modified z-score**, an
adopted robust outlier method [4]. For observed values \(x_i\), median
\(\tilde{x}\), and median absolute deviation (MAD):

\[
\operatorname{MAD}=\operatorname{median}\left(|x_i-\tilde{x}|\right)
\]

\[
M_i=\frac{0.6745(x_i-\tilde{x})}{\operatorname{MAD}}
\]

An observation is flagged as a potential outlier when:

\[
|M_i|>3.5
\]

These flags trigger manual review and never cause automatic exclusion. The
calculation is genus-wide and is therefore an anomaly screen, not a claim that
all *Bradyrhizobium* species have the same biological genome-size or GC
distribution.

Missing values are not treated as observations. If MAD is zero, non-median
values are considered unscored and require manual review; they must not be
silently treated as non-outliers.

## 6. Sequence-based quality control

### 6.1 CheckM2

CheckM2 completeness and contamination estimates are **adopted** as two
separate dimensions; no unsupported composite quality score will be invented.
CheckM2 was developed as a lineage-general machine-learning approach for
microbial genome quality estimation and may also be applied to isolate genomes
[5].

The decision thresholds are **adapted project criteria**:

- **Primary sequence-QC candidate:** completeness >=95% and contamination
  <=5%.
- **Secondary sequence-QC candidate:** completeness >=90% and contamination
  <=10%, but not satisfying both primary thresholds.
- **Exclude from the main pangenome:** completeness <90% or contamination
  >10%.

These rules are stricter than the MIMAG high-quality criterion in completeness
and use a similar contamination boundary. MIMAG was designed for MAGs and SAGs,
so it is used as a reference standard rather than an exact universal definition
for this mixed assembly collection [6].

The CheckM2 software version, compatible database version, database checksum,
command, input fingerprints, and raw `quality_report.tsv` must be preserved.

### 6.2 GUNC

GUNC is **adopted** to detect lineage-inconsistent chimerism using the genomic
distribution and taxonomic assignment of genes [7]. The native GUNC output,
including `pass.GUNC`, clade separation score, contamination portion, and
reference representation score, must be retained.

A GUNC failure, unscored result, or low-confidence reference representation
triggers manual review rather than automatic exclusion. GUNC cannot reliably
distinguish every horizontally transferred region from contamination. This is
especially important when evaluating symbiosis islands that may contain
legitimate horizontally transferred `nif`, `nod`, or `fix` genes [8].

A confirmed contamination decision must integrate GUNC, CheckM2, contig-level
taxonomy, assembly fragmentation, and genome context. The evidence and reviewer
reason must be recorded.

## 7. Two-stage taxonomy policy

Taxonomy is evaluated independently from sequence-quality estimates. NCBI and
GTDB assignments are both preserved; one must not silently overwrite the
other.

### 7.1 Stage 1: GTDB R220 screening

All available genomes are screened with GTDB R220 on a free high-memory
platform such as UseGalaxy.eu. R220 is a screening release, not the final
taxonomic authority for the project.

A genome enters targeted R232 review when one or more of the following apply:

- the GTDB genus is not `g__Bradyrhizobium`;
- the result is unclassified;
- GTDB-Tk reports warnings;
- no species assignment satisfies the release-specific ANI/AF criteria;
- NCBI and GTDB genus assignments disagree;
- metadata outlier flags are present;
- the assembly is strongly fragmented;
- another QC result makes the classification unreliable.

The project does not define a universal ANI boundary for genus membership.
Species decisions use GTDB's release-specific representative radius and
alignment-fraction policy [9].

### 7.2 Stage 2: GTDB R232 adjudication

The flagged subset is rerun against GTDB R232 on cloud or HPC infrastructure.
R232 contains substantially more genomes and species clusters than R220 and
uses an updated reference package and GTDB-Tk workflow [10-12].

R220 results must never be overwritten by R232 results. The minimum retained
fields are:

```text
gtdb_r220_classification
gtdb_r220_closest_reference
gtdb_r220_ani
gtdb_r220_af
gtdb_r220_warnings
gtdb_r232_classification
gtdb_r232_closest_reference
gtdb_r232_ani
gtdb_r232_af
gtdb_r232_warnings
taxonomy_release_concordance
taxonomy_review_decision
taxonomy_review_reason
```

An R220/R232 disagreement triggers manual review and does not automatically
exclude a genome. Only `taxonomy_review_decision`, supported by recorded
evidence, may change cohort membership.

For a publication-quality final analysis, the preferred endpoint is to run
R232 on the entire final cohort if resources permit. Until then, R220 provides
a consistent whole-cohort screen and R232 provides targeted adjudication.

## 8. Final cohort decisions

- **Primary:** available, metadata-resolved, satisfies primary contiguity and
  primary CheckM2 criteria, and has no unresolved GUNC or taxonomy anomaly.
- **Synteny candidate:** primary and also satisfies the synteny contiguity
  thresholds.
- **Secondary:** biologically usable but fragmented, satisfies only secondary
  CheckM2 criteria, or still requires a restricted downstream analysis.
- **Manual review:** metadata outlier, missing required metadata, GUNC
  failure/unscored result, taxonomy discordance, or another unresolved anomaly.
- **Exclude from main pangenome:** no current assembly member, CheckM2 below
  the exclusion boundary, confirmed contamination, or confirmed placement
  outside the target genus.
- **Pending technical:** download or tool failure. A technical failure must
  never be represented as biological exclusion.

All exclusion decisions require a machine-readable reason and retain the input
measurements that produced the decision.

## 9. Execution and resource policy

The user's laptop is suitable for metadata processing but not for the current
GTDB R232 bacterial workflow. Current GTDB-Tk documentation reports roughly
140 GB RAM and 100 GB storage for bacterial classification [13].

Execution therefore follows staged checkpoints:

1. Run a 10-genome pilot spanning normal, fragmented, metadata-outlier,
   broken-pair, and GCF-to-GCA fallback cases.
2. Record runtime, peak RAM, CPU allocation, tool/database versions, failures,
   and output schemas.
3. Run a 100-genome batch only after the pilot passes.
4. Run the full cohort in controlled batches only after resource and output
   validation.
5. Download and verify the compact TSV/JSON/log/manifest artifacts before
   deleting temporary cloud data.

UseGalaxy.eu is the preferred no-cost environment for the R220 whole-cohort
screen. Cloud or institutional HPC is reserved for R232 adjudication or a final
full-cohort R232 run. Platform-provided tool and database versions must be
recorded exactly.

## 10. Artifact contract

Planned tracked artifacts:

```text
docs/nif_pangenome_methodology_v1.md
config/sequence_qc_policy.json
data/metadata/bradyrhizobium_download_manifest.tsv
data/metadata/bradyrhizobium_sequence_qc.tsv
data/metadata/bradyrhizobium_manual_review.tsv
reports/task2_sequence_qc_summary.json
reports/task2_sequence_qc_run_manifest.json
```

Raw FASTA files, tool databases, and large raw tool directories are not stored
in Git. Their integrity and provenance are represented by checksums and run
manifests.

Every run manifest must record:

- exact accession-version input set and its SHA-256 fingerprint;
- software and dependency versions;
- reference database release and checksum;
- full command and relevant environment settings;
- start/end timestamps and execution platform;
- sample counts, success/failure counts, and retry history;
- output paths, output fingerprints, and schema versions.

Writers must be deterministic, use atomic replacement, and refuse silent
overwrite unless an explicit overwrite option is supplied.

## 11. Scope boundaries for later methodology versions

This version locks genome acquisition, metadata screening, sequence QC, and
taxonomy review. It does not yet lock:

- uniform genome annotation;
- gene-family clustering;
- formal core/accessory/rare frequency thresholds;
- pangenome openness modelling;
- `nif`, `nod`, and `fix` HMM panels and cutoffs;
- synteny distance/orientation rules;
- core-genome or symbiosis-gene phylogenies;
- horizontal gene transfer inference.

Each later stage requires its own audit, evidence-backed method design, artifact
contract, tests, and approval before execution on the research cohort.

## 12. Limitations

1. Metadata outliers may represent genuine biology or incomplete assemblies;
   they are not diagnoses.
2. CheckM2 estimates are model- and database-version dependent.
3. GUNC may confuse horizontally transferred regions with contamination and
   may be less informative when reference representation is poor.
4. GTDB taxonomy evolves between releases; species names, representatives,
   ANI matches, and placements may change.
5. R220 screening and targeted R232 adjudication are intentionally separate
   evidence layers. They must not be merged into an undocumented hybrid label.
6. Fragmentation can create false gene absence and disrupt apparent synteny.
7. Final biological conclusions about nitrogen fixation require gene-content,
   gene-context, and phylogenetic evidence; genome QC alone is insufficient.

## References

1. NCBI Datasets. Download large genome data packages.
   https://www.ncbi.nlm.nih.gov/datasets/docs/v2/how-tos/genomes/large-download/
2. Gurevich A, et al. QUAST: quality assessment tool for genome assemblies.
   *Bioinformatics*. 2013;29:1072-1075.
   https://doi.org/10.1093/bioinformatics/btt086
3. Tonkin-Hill G, et al. Producing polished prokaryotic pangenomes with the
   Panaroo pipeline. *Genome Biology*. 2020;21:180.
   https://doi.org/10.1186/s13059-020-02090-4
4. NIST/SEMATECH. Detection of Outliers: modified z-score criterion.
   https://itl.nist.gov/div898/handbook/eda/section3/eda35h.htm
5. Chklovski A, et al. CheckM2: a rapid, scalable and accurate tool for
   assessing microbial genome quality using machine learning. *Nature
   Methods*. 2023;20:1203-1212.
   https://doi.org/10.1038/s41592-023-01940-w
6. Bowers RM, et al. Minimum information about a single amplified genome and a
   metagenome-assembled genome of bacteria and archaea. *Nature
   Biotechnology*. 2017;35:725-731.
   https://doi.org/10.1038/nbt.3893
7. Orakov A, et al. GUNC: detection of chimerism and contamination in
   prokaryotic genomes. *Genome Biology*. 2021;22:178.
   https://doi.org/10.1186/s13059-021-02393-0
8. GUNC documentation. Output and limitations.
   https://grp-bork.embl-community.io/gunc/output.html
9. Chaumeil PA, et al. GTDB-Tk: a toolkit to classify genomes with the Genome
   Taxonomy Database. *Bioinformatics*. 2020;36:1925-1927.
   https://doi.org/10.1093/bioinformatics/btz848
10. GTDB R220 statistics.
    https://gtdb.ecogenomic.org/stats/r220
11. GTDB R232 statistics.
    https://gtdb.ecogenomic.org/stats/r232
12. GTDB-Tk change log.
    https://ecogenomics.github.io/GTDBTk/changelog.html
13. GTDB-Tk installation and hardware requirements.
    https://ecogenomics.github.io/GTDBTk/installing/index.html
