# HISAT-genotype


Please see the official website for HISATgenotype:
https://daehwankimlab.github.io/hisat-genotype/

*This directory no longer contains HISAT2. Please see the link below if you are looking for the newest HISAT2 release*
https://daehwankimlab.github.io/hisat2/

## Current Release Version
v1.3.2 - Update to patch critical extract reads error

## Previous Releases
v1.3.1 - Update database locations and handling, and general stability fixes - BROKEN

v1.3.0 - Python 3 version

## Overview
HISAT-genotype is a next-generation genomic analysis software platform capable of assembling and genotyping human genes and genomic regions. Thie software leverages HISAT2s graph FM index and graph alignemnt algorithm to align reads to a specially constructed graph genome. An Expectation-Maximization (EM) algorithm finds the maximum likelihood estimates for each gene allele and a guided de Bruijn graph is used to construct the allele sequences.

## HLA Pipeline Performance Benchmark (RedStation, 2026-06-03)

5 whole-genome (30x WGS, ~30-37 GB BAM) samples run in parallel through the full
HLA pipeline (HLA-LA -> HISAT-genotype -> HLA-FAST -> OptiType -> consensus
consolidation). Reference data mounted read-only; outputs written to job workdir.
Sample identities anonymized.

| Genome | BAM size | Node placement | Wall-clock | CPU used | Avg cores | Peak mem* |
|--------|----------|----------------|------------|----------|-----------|-----------|
| 1      | 30 GB    | dedicated node | 72 min     | 24.5 CPU-h | 20.4    | 60 GB     |
| 2      | 36 GB    | dedicated node | 78 min     | 26.0 CPU-h | 20.0    | 69 GB     |
| 3      | 32 GB    | shared (3/node)| 150 min    | 30.3 CPU-h | 12.1    | 62 GB     |
| 4      | 37 GB    | shared (3/node)| 168 min    | 31.7 CPU-h | 11.3    | 69 GB     |
| 5      | 31 GB    | shared (3/node)| 183 min    | 26.8 CPU-h |  8.8    | 75 GB     |

\* Peak memory is cgroup `memory.current`; roughly half is reclaimable BAM page
cache. True resident set (dominated by the HLA-LA PRG MHC graph) is ~28-30 GB.

### Findings
- **Compute per genome is ~constant (~25-32 CPU-hours)** regardless of placement.
- **Wall-clock is set by available cores.** A dedicated node (~20 cores) finishes
  in ~75 min; three pods sharing one node throttle to ~9-12 cores each and take
  150-183 min (~2-2.5x penalty) -- contention is worst during the HLA-LA phase.
- Pipeline is **CPU-bound**, not I/O (BAM download was ~7 min of total).
- HLA-LA reads the full BAM only to extract MHC-region reads (~234 MB paired
  FASTQ); typing runs on those.

### Resource recommendation
- Run **1 HLA pod per node** (pod anti-affinity) or cap concurrency so requested
  cores <= node cores. Budget **~20 cores and ~30 GB RAM** per pod.
- Throughput scales with total cores, not pod count.
