#!/usr/bin/env python3
"""
Regenerate ALL _nuc.msf files from _gen.msf for IMGT/HLA v3.61
to fix exon coordinate mismatches
"""

# Import the regeneration function
import sys
sys.path.insert(0, '/root/hisat-genotype')
from regenerate_nuc_msf import regenerate_nuc_msf
import os

# Configuration for v3.61
hla_db_dir = "/root/hisat-genotype/hisatgenotype_db/HLA"
hla_dat = f"{hla_db_dir}/hla.dat"
msf_dir = f"{hla_db_dir}/msf"

# ALL genes with exon coordinate mismatches in v3.61
genes_with_issues = [
    'A', 'B', 'C', 'DQB1', 'DRB1', 'DPA1', 'DPB1', 
    'DRB5', 'E', 'J', 'DRB3', 'DRB4'
]

# Reference alleles (from FASTA files)
reference_alleles = {
    'A': 'A*01:01:01:01',
    'B': 'B*07:02:01:01',
    'C': 'C*01:02:01:01',
    'DQB1': 'DQB1*05:01:01:01',
    'DRB1': 'DRB1*01:01:01:01',
    'DPA1': 'DPA1*01:03:01:01',
    'DPB1': 'DPB1*01:01:01:01',
    'DRB5': 'DRB5*01:01:01:01',
    'E': 'E*01:01:01:01',
    'J': 'J*01:01:01:01',
    'DRB3': 'DRB3*01:01:02:01',
    'DRB4': 'DRB4*01:01:01:01',
}

print("="*70)
print("REGENERATING _nuc.msf FILES FOR v3.61")
print("="*70)
print(f"Database: {hla_db_dir}")
print(f"Processing {len(genes_with_issues)} genes with exon mismatches")
print()

success_count = 0
failed_genes = []

for gene in genes_with_issues:
    gen_msf = f"{msf_dir}/{gene}_gen.msf"
    nuc_msf = f"{msf_dir}/{gene}_nuc.msf"
    nuc_backup = f"{msf_dir}/{gene}_nuc.msf.v361_original"
    ref_allele = reference_alleles.get(gene)
    
    if not ref_allele:
        print(f"✗ {gene}: No reference allele defined")
        failed_genes.append(gene)
        continue
    
    if not os.path.exists(gen_msf):
        print(f"✗ {gene}: {gen_msf} not found")
        failed_genes.append(gene)
        continue
    
    # Backup original _nuc.msf
    if os.path.exists(nuc_msf) and not os.path.exists(nuc_backup):
        os.rename(nuc_msf, nuc_backup)
        print(f"  Backed up original to {os.path.basename(nuc_backup)}")
    
    try:
        if regenerate_nuc_msf(gene, gen_msf, hla_dat, nuc_msf, ref_allele):
            success_count += 1
            print(f"✓ {gene}: Successfully regenerated")
        else:
            failed_genes.append(gene)
            print(f"✗ {gene}: Regeneration failed")
    except Exception as e:
        failed_genes.append(gene)
        print(f"✗ {gene}: ERROR - {e}")

print()
print("="*70)
print(f"RESULTS: {success_count}/{len(genes_with_issues)} genes successful")
if failed_genes:
    print(f"Failed genes: {', '.join(failed_genes)}")
print("="*70)
