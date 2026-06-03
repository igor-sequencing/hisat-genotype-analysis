#!/usr/bin/env python3
"""
Split DRB345 MSF files into individual DRB3, DRB4, DRB5 files.
"""

import sys
import os
from collections import OrderedDict

def split_msf(input_file, output_dir):
    """Split combined DRB345 MSF file into individual gene files."""

    print(f"Reading {input_file}...")

    with open(input_file, 'r') as f:
        lines = f.readlines()

    # Find header end
    header_end = None
    for i, line in enumerate(lines):
        if line.strip() == '//':
            header_end = i
            break

    if header_end is None:
        print("Error: Could not find '//' marker")
        return False

    # Parse header to get allele info
    alleles_by_gene = {'DRB3': OrderedDict(), 'DRB4': OrderedDict(), 'DRB5': OrderedDict()}
    header_lines = lines[:header_end]

    msf_line = header_lines[2]  # Line 3 has MSF info

    for line in header_lines:
        if line.strip().startswith('Name:'):
            parts = line.strip().split()
            allele_name = parts[1]
            gene = allele_name.split('*')[0]

            if gene in alleles_by_gene:
                alleles_by_gene[gene][allele_name] = {
                    'len': int(parts[3]),
                    'check': int(parts[5]),
                    'weight': float(parts[7])
                }

    # Count alleles per gene
    for gene, alleles in alleles_by_gene.items():
        print(f"  {gene}: {len(alleles)} alleles")

    # Parse alignment section
    alignment_lines = lines[header_end:]

    # Group alignment lines by allele
    allele_sequences = {}
    for line in alignment_lines:
        line = line.strip()
        if not line or line == '//':
            continue
        parts = line.split()
        if len(parts) >= 2:
            allele_name = parts[0]
            seq_part = ''.join(parts[1:])
            if allele_name not in allele_sequences:
                allele_sequences[allele_name] = []
            allele_sequences[allele_name].append(seq_part)

    # Write separate files for each gene
    base_name = os.path.basename(input_file)
    file_type = base_name.replace('DRB345', '').replace('.msf', '')  # _nuc or _prot

    for gene in ['DRB3', 'DRB4', 'DRB5']:
        if not alleles_by_gene[gene]:
            print(f"  Skipping {gene} (no alleles found)")
            continue

        output_file = os.path.join(output_dir, f"{gene}{file_type}.msf")
        print(f"  Writing {output_file}...")

        with open(output_file, 'w') as f:
            # Write header
            f.write("!!NA_MULTIPLE_ALIGNMENT\n\n")
            f.write(msf_line)
            f.write("\n")

            # Write allele info
            for allele_name, info in alleles_by_gene[gene].items():
                f.write(f" Name: {allele_name:20s} Len: {info['len']:5d}  Check: {info['check']:5d}  Weight:  {info['weight']:.2f}\n")

            f.write("\n//\n\n")

            # Write alignment
            for allele_name in alleles_by_gene[gene].keys():
                if allele_name in allele_sequences:
                    seq_blocks = allele_sequences[allele_name]
                    # Write sequence in blocks of 50 characters per line
                    full_seq = ''.join(seq_blocks)
                    for i in range(0, len(full_seq), 50):
                        block = full_seq[i:i+50]
                        f.write(f"{allele_name:20s} {block}\n")
                    f.write("\n")

    return True

def main():
    msf_dir = "indicies/hisatgenotype_db/HLA/msf"

    # Split both _nuc and _prot files
    for file_type in ['_nuc', '_prot']:
        input_file = os.path.join(msf_dir, f"DRB345{file_type}.msf")

        if os.path.exists(input_file):
            print(f"\nProcessing {input_file}...")
            success = split_msf(input_file, msf_dir)
            if success:
                print(f"  ✓ Successfully split DRB345{file_type}.msf")
            else:
                print(f"  ✗ Failed to split DRB345{file_type}.msf")
        else:
            print(f"Warning: {input_file} not found")

    print("\nDone!")

if __name__ == "__main__":
    main()
