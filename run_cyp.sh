#!/bin/bash

samtools view -b -L /mnt/data/refData/hla-tools/hisat-genotype/cyp.bed $1 >cyp.bam
samtools index cyp.bam

samtools sort -n -o - cyp.bam | samtools fastq -1 cyp_R1.fq.gz -2 cyp_R2.fq.gz -s cyp_s.fq.gz -0 cyp_0.fq.gz -

hisatgenotype \
  -x genotype_genome \
  -z /mnt/data/refData/hla-tools/hisat-genotype/indicies \
  --base CYP \
  -1 cyp_R1.fq.gz \
  -2 cyp_R2.fq.gz \
  --out-dir hisat.output.2 
