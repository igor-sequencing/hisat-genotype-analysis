#!/usr/bin/env python3
"""
Validate MSF database files for all HLA genes
Checks for file existence, reference allele consistency, and exon coordinate matching
"""

import sys
import os
import re

# Add the hisatgenotype modules to the path
sys.path.insert(0, '/root/hisat-genotype')
sys.path.insert(0, '/root/hisat-genotype/hisatgenotype_modules')

import hisatgenotype_typing_common as typing_common

# List of all 32 HLA genes
HLA_GENES = [
    "A", "B", "C", "DQA1", "DQB1", "DRB1",
    "DMA", "DMB", "DOA", "DOB", "DPA1", "DPB1", "DRA", "DRB5",
    "E", "F", "G", "H", "J", "K", "L",
    "MICA", "MICB", "P", "TAP1", "TAP2", "V", "W", "Y",
    "DRB3", "DRB4", "ClassI"
]

# Genes that have spliced exons (should have _nuc.msf)
SPLICED_GENES = [
    "A", "B", "C", "DQA1", "DQB1", "DRB1",
    "DMA", "DMB", "DOA", "DOB", "DPA1", "DPB1", "DRA", "DRB5",
    "E", "F", "G", "H", "J", "K", "L",
    "MICA", "MICB", "TAP1", "TAP2",
    "DRB3", "DRB4", "ClassI"
]

def read_fasta_file(fname):
    """Read a FASTA file and return the reference allele name"""
    if not os.path.exists(fname):
        return None

    with open(fname, 'r') as f:
        for line in f:
            if line.startswith('>'):
                # Extract allele name from header
                header = line[1:].strip()
                # Extract gene*allele pattern
                match = re.search(r'(HLA-)?([A-Z0-9]+)\*([0-9:]+)', header)
                if match:
                    gene = match.group(2)
                    allele = match.group(3)
                    return f"{gene}*{allele}"
                return header
    return None

def find_reference_in_msf(msf_fname, ref_allele_from_fasta):
    """Check if the reference allele from FASTA exists in MSF file"""
    if not os.path.exists(msf_fname):
        return None, []

    try:
        # read_MSF_file requires full_alleles dict - pass empty dict for validation
        full_alleles = {}
        names, seqs = typing_common.read_MSF_file(msf_fname, full_alleles, "", "")

        # Check for exact match
        if ref_allele_from_fasta in names:
            return True, list(names.keys())

        # Check for partial match (different version suffixes)
        base_ref = ref_allele_from_fasta.split('.')[0]  # Remove version suffix
        for name in names.keys():
            base_name = name.split('.')[0]
            if base_name == base_ref or base_name == ref_allele_from_fasta or name == ref_allele_from_fasta:
                return f"PARTIAL:{name}", list(names.keys())

        # Check for prefix match (e.g., "MICA*008:04" matches "MICA*008:04:01")
        for name in names.keys():
            if ref_allele_from_fasta.startswith(name + ":") or name.startswith(ref_allele_from_fasta + ":"):
                return f"PARTIAL:{name}", list(names.keys())

        return False, list(names.keys())
    except Exception as e:
        return f"ERROR:{str(e)}", []

def extract_exon_positions(gene, ref_allele_name):
    """Extract exon positions from the .dat file"""
    dat_file = "hisatgenotype_db/HLA/hla.dat"

    if not os.path.exists(dat_file):
        return None

    exons = []
    current_allele = None
    look_exon_num = False

    # Read the .dat file to find exon positions for the reference allele
    for line in open(dat_file):
        if line.startswith("DE"):
            # Extract allele name from DE line
            if not line.split()[1][-1].isdigit():
                allele_name = line.split()[1][:-1]
            else:
                allele_name = line.split()[1]

            # Remove HLA- prefix if present
            if allele_name.startswith("HLA-"):
                allele_name = allele_name[4:]

            current_allele = allele_name

        if not line.startswith("FT"):
            continue

        if line.find("exon") != -1:
            look_exon_num = True
            # Check if this is the reference allele (or a prefix match)
            if current_allele and (current_allele == ref_allele_name or
                                  current_allele.startswith(ref_allele_name + ":") or
                                  ref_allele_name.startswith(current_allele + ":")):
                exon_range = line.split()[2].split("..")
                exon_left = int(exon_range[0]) - 1
                exon_right = int(exon_range[1]) - 1
                exons.append([exon_left, exon_right])
        elif look_exon_num:
            look_exon_num = False

    return exons if exons else None

def create_map(seq):
    """Create a map from base pair position to MSF position"""
    seq_map = {}
    count = 0
    for i in range(len(seq)):
        bp = seq[i]
        if bp in '.EN~':
            continue
        seq_map[count] = i
        count += 1
    return seq_map

def compare_exon_coordinates(gene, gen_fname, nuc_fname, ref_allele_name):
    """Compare exon coordinates between genomic and nucleotide MSF files"""
    if not os.path.exists(gen_fname) or not os.path.exists(nuc_fname):
        return None

    try:
        # Read both files
        full_alleles = {}
        gen_names, gen_seqs = typing_common.read_MSF_file(gen_fname, full_alleles, "", "")
        nuc_names, nuc_seqs = typing_common.read_MSF_file(nuc_fname, full_alleles)

        # Find the reference allele in both files
        ref_gen = None
        ref_nuc = None

        # Try exact match first
        if ref_allele_name in gen_names:
            ref_gen = ref_allele_name
        if ref_allele_name in nuc_names:
            ref_nuc = ref_allele_name

        # Try prefix matching if exact match not found
        if not ref_gen:
            for name in gen_names.keys():
                if name.startswith(ref_allele_name + ":") or ref_allele_name.startswith(name + ":"):
                    ref_gen = name
                    break

        if not ref_nuc:
            for name in nuc_names.keys():
                if name.startswith(ref_allele_name + ":") or ref_allele_name.startswith(name + ":"):
                    ref_nuc = name
                    break

        if not ref_gen or not ref_nuc:
            return "NO_COMMON_ALLELES"

        # Get sequences
        gen_seq = gen_seqs[gen_names[ref_gen]]
        nuc_seq = nuc_seqs[nuc_names[ref_nuc]]

        # Get exon positions from .dat file
        exons = extract_exon_positions(gene, ref_allele_name)
        if exons is None or len(exons) == 0:
            return "NO_EXON_DATA"

        # Create position maps
        gen_seq_map = create_map(gen_seq)
        nuc_seq_map = create_map(nuc_seq)

        # Compare each exon
        exon_mismatches = []
        nuc_pos = 0
        for exon_i, (left, right) in enumerate(exons):
            # Get exon from genomic sequence
            if left not in gen_seq_map or right not in gen_seq_map:
                exon_mismatches.append(f"exon{exon_i}:OUT_OF_BOUNDS")
                continue

            gen_exon_start = gen_seq_map[left]
            gen_exon_end = gen_seq_map[right]
            gen_exon_len = gen_exon_end - gen_exon_start + 1

            # Get corresponding region from nucleotide sequence
            nuc_exon_end = nuc_pos + (right - left)
            if nuc_exon_end >= len(nuc_seq_map):
                exon_mismatches.append(f"exon{exon_i}:NUC_TOO_SHORT")
                break

            nuc_exon_start = nuc_seq_map[nuc_pos] if nuc_pos in nuc_seq_map else -1
            nuc_exon_end_pos = nuc_seq_map[nuc_exon_end] if nuc_exon_end in nuc_seq_map else -1

            if nuc_exon_start < 0 or nuc_exon_end_pos < 0:
                exon_mismatches.append(f"exon{exon_i}:NUC_MAP_ERROR")
                nuc_pos = nuc_exon_end + 1
                continue

            nuc_exon_len = nuc_exon_end_pos - nuc_exon_start + 1

            if gen_exon_len != nuc_exon_len:
                exon_mismatches.append(f"exon{exon_i}:gen={gen_exon_len},nuc={nuc_exon_len}")

            nuc_pos = nuc_exon_end + 1

        if exon_mismatches:
            return f"MISMATCH:{';'.join(exon_mismatches)}"
        else:
            return "MATCH"

    except Exception as e:
        return f"ERROR:{str(e)}"

def validate_all_genes():
    """Validate MSF files for all HLA genes"""

    print("=" * 100)
    print("HLA MSF FILE VALIDATION REPORT")
    print("=" * 100)
    print()

    base_dir = "hisatgenotype_db/HLA"

    results = {
        'genes_checked': 0,
        'all_files_present': 0,
        'missing_gen_msf': [],
        'missing_gen_fasta': [],
        'missing_nuc_msf': [],
        'missing_prot_msf': [],
        'ref_allele_mismatch': [],
        'exon_coord_mismatch': [],
        'errors': []
    }

    for gene in HLA_GENES:
        results['genes_checked'] += 1

        print(f"\n{'=' * 80}")
        print(f"Gene: {gene}")
        print(f"{'=' * 80}")

        # File paths (files use uppercase names)
        gen_msf = f"{base_dir}/msf/{gene}_gen.msf"
        gen_fasta = f"{base_dir}/fasta/{gene}_gen.fasta"
        nuc_msf = f"{base_dir}/msf/{gene}_nuc.msf"
        prot_msf = f"{base_dir}/msf/{gene}_prot.msf"

        # Check file existence
        gen_msf_exists = os.path.exists(gen_msf)
        gen_fasta_exists = os.path.exists(gen_fasta)
        nuc_msf_exists = os.path.exists(nuc_msf)
        prot_msf_exists = os.path.exists(prot_msf)

        print(f"\nFile Existence:")
        print(f"  {gene}_gen.msf:    {'✓ EXISTS' if gen_msf_exists else '✗ MISSING'}")
        print(f"  {gene}_gen.fasta: {'✓ EXISTS' if gen_fasta_exists else '✗ MISSING'}")
        print(f"  {gene}_nuc.msf:   {'✓ EXISTS' if nuc_msf_exists else '✗ MISSING (expected)' if gene not in SPLICED_GENES else '✗ MISSING'}")
        print(f"  {gene}_prot.msf:  {'✓ EXISTS' if prot_msf_exists else '✗ MISSING'}")

        # Track missing files
        if not gen_msf_exists:
            results['missing_gen_msf'].append(gene)
        if not gen_fasta_exists:
            results['missing_gen_fasta'].append(gene)
        if gene in SPLICED_GENES and not nuc_msf_exists:
            results['missing_nuc_msf'].append(gene)
        if not prot_msf_exists:
            results['missing_prot_msf'].append(gene)

        # Check reference allele consistency
        if gen_fasta_exists and gen_msf_exists:
            ref_from_fasta = read_fasta_file(gen_fasta)
            ref_in_msf, msf_alleles = find_reference_in_msf(gen_msf, ref_from_fasta)

            print(f"\nReference Allele Consistency:")
            print(f"  FASTA reference: {ref_from_fasta}")
            print(f"  MSF alleles: {len(msf_alleles)} total")

            if ref_in_msf == True:
                print(f"  ✓ Reference allele found in MSF (exact match)")
            elif isinstance(ref_in_msf, str) and ref_in_msf.startswith("PARTIAL:"):
                matched_name = ref_in_msf.split(":", 1)[1]
                print(f"  ⚠ Reference allele found with different suffix: {matched_name}")
                results['ref_allele_mismatch'].append(f"{gene}:VERSION_MISMATCH")
            elif ref_in_msf == False:
                print(f"  ✗ Reference allele NOT found in MSF")
                print(f"     First few MSF alleles: {', '.join(msf_alleles[:5])}")
                results['ref_allele_mismatch'].append(f"{gene}:NOT_FOUND")
            else:
                print(f"  ✗ ERROR: {ref_in_msf}")
                results['errors'].append(f"{gene}:ref_check:{ref_in_msf}")

        # Check exon coordinate consistency
        if gene in SPLICED_GENES and gen_msf_exists and nuc_msf_exists and gen_fasta_exists:
            ref_from_fasta = read_fasta_file(gen_fasta)
            exon_check = compare_exon_coordinates(gene, gen_msf, nuc_msf, ref_from_fasta)

            print(f"\nExon Coordinate Consistency:")
            if exon_check == "MATCH":
                print(f"  ✓ Exon coordinates match between _gen.msf and _nuc.msf")
            elif exon_check and exon_check.startswith("MISMATCH"):
                details = exon_check.split(":", 1)[1]
                print(f"  ✗ Exon coordinate MISMATCH: {details}")
                results['exon_coord_mismatch'].append(f"{gene}:{details}")
            elif exon_check == "NO_COMMON_ALLELES":
                print(f"  ✗ No common alleles found between _gen.msf and _nuc.msf")
                results['exon_coord_mismatch'].append(f"{gene}:NO_COMMON_ALLELES")
            elif exon_check == "NO_EXON_DATA":
                print(f"  ⚠ No exon data found in .dat file")
            elif exon_check and exon_check.startswith("ERROR"):
                error_msg = exon_check.split(":", 1)[1]
                print(f"  ✗ ERROR: {error_msg}")
                results['errors'].append(f"{gene}:exon_check:{error_msg}")

        # Count genes with all required files
        if gen_msf_exists and gen_fasta_exists:
            if gene in SPLICED_GENES:
                if nuc_msf_exists:
                    results['all_files_present'] += 1
            else:
                results['all_files_present'] += 1

    # Print summary
    print("\n\n" + "=" * 100)
    print("VALIDATION SUMMARY")
    print("=" * 100)
    print(f"\nGenes checked: {results['genes_checked']}")
    print(f"Genes with all required files: {results['all_files_present']}")

    if results['missing_gen_msf']:
        print(f"\n✗ Missing _gen.msf files ({len(results['missing_gen_msf'])}):")
        for gene in results['missing_gen_msf']:
            print(f"    - {gene}")

    if results['missing_gen_fasta']:
        print(f"\n✗ Missing _gen.fasta files ({len(results['missing_gen_fasta'])}):")
        for gene in results['missing_gen_fasta']:
            print(f"    - {gene}")

    if results['missing_nuc_msf']:
        print(f"\n✗ Missing _nuc.msf files for spliced genes ({len(results['missing_nuc_msf'])}):")
        for gene in results['missing_nuc_msf']:
            print(f"    - {gene}")

    if results['missing_prot_msf']:
        print(f"\n✗ Missing _prot.msf files ({len(results['missing_prot_msf'])}):")
        for gene in results['missing_prot_msf']:
            print(f"    - {gene}")

    if results['ref_allele_mismatch']:
        print(f"\n⚠ Reference allele mismatches ({len(results['ref_allele_mismatch'])}):")
        for issue in results['ref_allele_mismatch']:
            print(f"    - {issue}")

    if results['exon_coord_mismatch']:
        print(f"\n✗ Exon coordinate mismatches ({len(results['exon_coord_mismatch'])}):")
        for issue in results['exon_coord_mismatch']:
            print(f"    - {issue}")

    if results['errors']:
        print(f"\n✗ Errors during validation ({len(results['errors'])}):")
        for error in results['errors']:
            print(f"    - {error}")

    print("\n" + "=" * 100)

    return results

if __name__ == "__main__":
    validate_all_genes()
