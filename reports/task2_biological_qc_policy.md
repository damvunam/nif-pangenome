Task 2 — Biological QC policy

Scope

This policy separates metadata screening from later sequence-based quality
assessment. Metadata screening must not be described as final genome QC.

Cohorts

Primary pangenome candidate: a current assembly with no unresolved
metadata anomaly, no more than 500 contigs, and contig N50 of at least 20 kb.

Synteny candidate: a primary candidate with no more than 200 contigs and
contig N50 of at least 50 kb.

Secondary candidate: an available assembly that is too fragmented for the
primary analysis or requires manual metadata review.

Unavailable: no current GCA or GCF member remains for the physical
assembly.

The contiguity thresholds are project-specific safeguards, not universal
definitions of genome quality. N50 is never used as a stand-alone exclusion
criterion because it can increase despite misassembly (Gurevich et al.
2013). Fragmentation and
annotation errors can distort core and accessory gene estimates
(Tonkin-Hill et al. 2020).

Representative selection

For each canonical physical assembly:

Prefer a current RefSeq (GCF) member.

If the GCF member is suppressed but its GenBank (GCA) partner is current,
use the GCA record and retain the replacement reason.

Keep a current broken-pair GCA as a candidate with an explicit traceability
flag.

Mark a physical assembly unavailable only if neither member is current.

NCBI notes that prokaryotic RefSeq assemblies may be suppressed because of
sequence, annotation, or metadata issues, while the paired GenBank assembly
may remain available (NCBI assembly status
documentation).

Later sequence QC

Metadata screening is followed by CheckM2 completeness and contamination
estimation. The primary cohort requires at least 95% completeness and no more
than 5% contamination. Assemblies with 90–<95% completeness or >5–10%
contamination remain secondary candidates; lower completeness or higher
contamination excludes them from the main pangenome. These project thresholds
are stricter than the >90% completeness and <5% contamination MIMAG
high-quality-draft benchmarks (Bowers et al.
2017). CheckM2 provides lineage-general
microbial genome quality estimates (Chklovski et al.
2023).

GUNC failures require contig-level review rather than automatic exclusion.
GUNC detects chimerism missed by marker-based methods (Orakov et al.
2021), but this project must avoid
confusing genuine horizontally transferred symbiosis regions with contamination.

Outliers and taxonomy

Genome size and GC percentage are flagged with a modified z-score threshold of
3.5. They do not cause automatic exclusion. Taxonomy must later be checked with
GTDB-Tk and ANI/manual verification before any genome is declared outside the
target genus.