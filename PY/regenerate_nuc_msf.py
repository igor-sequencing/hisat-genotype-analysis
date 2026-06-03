#!/usr/bin/env python3
"""
Regenerate _nuc.msf files from _gen.msf using hla.dat exon coordinates.

This fixes the gap pattern inconsistency between genomic and coding MSF alignments
by extracting exonic regions from the genomic MSF and preserving its gap structure.
"""

import sys
import os
import re
from collections import OrderedDict

def read_msf_file(msf_file):
    """Read MSF file and return header, allele names, and sequences."""
    with open(msf_file, 'r') as f:
        lines = f.readlines()

    # Find alignment start
    alignment_start = None
    for i, line in enumerate(lines):
        if line.strip() == '//':
            alignment_start = i + 1
            break

    if alignment_start is None:
        raise ValueError(f"No '//' marker found in {msf_file}")

    header = lines[:alignment_start]
    alignment_lines = lines[alignment_start:]

    # Extract allele names from header
    allele_info = OrderedDict()  # Preserve order
    for line in header:
        if line.strip().startswith('Name:'):
            parts = line.strip().split()
            allele_name = parts[1]
            length = int(parts[3])
            checksum = int(parts[5])
            allele_info[allele_name] = {'len': length, 'check': checksum}

    # Extract sequences
    sequences = {name: [] for name in allele_info.keys()}
    for line in alignment_lines:
        line = line.strip()
        if not line or line.startswith('//'):
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        allele_name = parts[0]
        if allele_name in sequences:
            # Join all sequence blocks on this line
            seq_parts = ''.join(parts[1:])
            sequences[allele_name].append(seq_parts)

    # Concatenate sequence blocks
    for name in sequences:
        sequences[name] = ''.join(sequences[name])

    return header, allele_info, sequences


def create_position_map(seq_with_gaps):
    """Create mapping from biological position (no gaps) to alignment position (with gaps)."""
    bio_to_aln = {}
    bio_pos = 0

    for aln_pos, char in enumerate(seq_with_gaps):
        if char not in '.~EN':  # Not a gap or special character
            bio_to_aln[bio_pos] = aln_pos
            bio_pos += 1

    return bio_to_aln


def extract_exon_coordinates(hla_dat, gene, ref_allele):
    """Extract exon coordinates for a specific gene from hla.dat."""
    exons = []
    in_target_allele = False

    with open(hla_dat, 'r') as f:
        for line in f:
            if line.startswith("DE"):
                # Check if this is our reference allele
                if not line.split()[1][-1].isdigit():
                    allele_name = line.split()[1][:-1]
                else:
                    allele_name = line.split()[1]

                # Remove HLA- prefix if present
                if allele_name.startswith("HLA-"):
                    allele_name = allele_name[4:]

                if allele_name == ref_allele:
                    in_target_allele = True
                else:
                    in_target_allele = False

            if not in_target_allele:
                continue

            if line.startswith("FT") and "exon" in line and ".." in line:
                # Extract exon range
                parts = line.split()
                for part in parts:
                    if ".." in part:
                        exon_range = part.split("..")
                        exon_left = int(exon_range[0]) - 1  # Convert to 0-indexed
                        exon_right = int(exon_range[1]) - 1
                        exons.append([exon_left, exon_right])
                        break

    return exons


def align_partial_to_pattern(partial_seq_nogap, reference_with_gaps):
    """
    Align a partial allele (no gaps) to match the gap pattern of a reference.
    Insert gaps at the same positions as in the reference.
    """
    result = []
    bio_pos = 0

    for char in reference_with_gaps:
        if char in '.~':
            # Insert gap
            result.append(char)
        else:
            # Insert next biological base
            if bio_pos < len(partial_seq_nogap):
                result.append(partial_seq_nogap[bio_pos])
                bio_pos += 1
            else:
                # Ran out of sequence, pad with gaps
                result.append('.')

    # If we have leftover sequence, append it
    while bio_pos < len(partial_seq_nogap):
        result.append(partial_seq_nogap[bio_pos])
        bio_pos += 1

    return ''.join(result)


def regenerate_nuc_msf(gene, gen_msf_path, hla_dat_path, output_path, ref_allele):
    """
    Regenerate _nuc.msf from _gen.msf using exon coordinates from hla.dat.
    Also includes partial alleles from original _nuc.msf, converted to new gap pattern.
    """
    print(f"\n=== Processing {gene} ===")

    # Read genomic MSF
    print(f"Reading {gen_msf_path}...")
    header, allele_info, gen_sequences = read_msf_file(gen_msf_path)

    # Get exon coordinates
    print(f"Extracting exon coordinates for {ref_allele}...")
    exons = extract_exon_coordinates(hla_dat_path, gene, ref_allele)
    if not exons:
        print(f"ERROR: No exons found for {gene} / {ref_allele}")
        return False

    print(f"Found {len(exons)} exons: {exons[:3]}...")  # Show first 3

    # Get reference sequence and create position map
    ref_seq = gen_sequences[ref_allele]
    bio_to_aln = create_position_map(ref_seq)

    print(f"Reference sequence: {len(ref_seq)} alignment chars, {len(bio_to_aln)} biological bp")

    # Extract exonic regions for all alleles from genomic MSF
    nuc_sequences = {}
    nuc_lengths = []

    for allele_name, gen_seq in gen_sequences.items():
        if len(gen_seq) != len(ref_seq):
            print(f"  Warning: {allele_name} has different length ({len(gen_seq)} vs {len(ref_seq)}), skipping")
            continue

        # Extract and concatenate exons
        exonic_seq = []
        for exon_left, exon_right in exons:
            # Map biological coordinates to alignment positions
            if exon_left not in bio_to_aln or exon_right not in bio_to_aln:
                print(f"  Warning: {allele_name} missing exon coordinates, skipping")
                break
            aln_left = bio_to_aln[exon_left]
            aln_right = bio_to_aln[exon_right]

            # Extract exon with gaps from genomic sequence
            exon_seq = gen_seq[aln_left:aln_right + 1]
            exonic_seq.append(exon_seq)
        else:
            # All exons extracted successfully
            nuc_seq = ''.join(exonic_seq)
            nuc_sequences[allele_name] = nuc_seq
            nuc_lengths.append(len(nuc_seq))

    if not nuc_sequences:
        print(f"ERROR: No sequences extracted for {gene}")
        return False

    # Check all sequences have same length
    if len(set(nuc_lengths)) != 1:
        print(f"ERROR: Sequences have different lengths: {set(nuc_lengths)}")
        return False

    nuc_length = nuc_lengths[0]
    print(f"Extracted {len(nuc_sequences)} full alleles with length {nuc_length}")

    # Now add partial alleles from original _nuc.msf
    original_nuc = gen_msf_path.replace('_gen.msf', '_nuc.msf')
    # Check for .original backup first (in case we already replaced the file)
    if os.path.exists(original_nuc + '.original'):
        original_nuc = original_nuc + '.original'
    if os.path.exists(original_nuc):
        print(f"\nAdding partial alleles from {os.path.basename(original_nuc)}...")
        _, orig_allele_info, orig_sequences = read_msf_file(original_nuc)

        # Find partial alleles (those NOT in genomic MSF)
        partial_alleles = set(orig_sequences.keys()) - set(gen_sequences.keys())
        print(f"Found {len(partial_alleles)} partial alleles")

        # Get reference allele from regenerated sequences to use as gap pattern template
        ref_regenerated = nuc_sequences[ref_allele]

        added_count = 0
        for allele_name in partial_alleles:
            orig_seq = orig_sequences[allele_name]
            # Remove gaps to get biological sequence
            bio_seq = orig_seq.replace('.', '').replace('~', '').replace('E', '').replace('N', '')

            # Align to new gap pattern (may be longer if allele has insertions)
            new_seq = align_partial_to_pattern(bio_seq, ref_regenerated)
            nuc_sequences[allele_name] = new_seq
            added_count += 1

        print(f"Successfully added {added_count} partial alleles")

        # Find maximum length (some alleles may have insertions)
        all_lengths = [len(seq) for seq in nuc_sequences.values()]
        max_length = max(all_lengths)

        if max_length > nuc_length:
            print(f"  Note: Some alleles have insertions, max length is {max_length} (reference: {nuc_length})")
            print(f"  Padding all sequences to {max_length} bp with gaps")
            # Pad all shorter sequences with gaps
            for allele_name in nuc_sequences:
                seq = nuc_sequences[allele_name]
                if len(seq) < max_length:
                    nuc_sequences[allele_name] = seq + '.' * (max_length - len(seq))
            nuc_length = max_length

        print(f"Total alleles: {len(nuc_sequences)} ({len(nuc_sequences) - added_count} full + {added_count} partial)")

    # Write new MSF file
    print(f"Writing {output_path}...")
    with open(output_path, 'w') as f:
        # Write MSF header
        f.write("!!NA_MULTIPLE_ALIGNMENT\n\n")
        f.write(f"   MSF: {nuc_length}  Type: N  Regenerated  Check: 0 ..\n\n")

        # Write Name entries
        for allele_name in nuc_sequences.keys():
            # Calculate simple checksum (sum of ord values mod 10000)
            checksum = sum(ord(c) for c in nuc_sequences[allele_name]) % 10000
            f.write(f" Name: {allele_name:<15} Len: {nuc_length:5}  Check: {checksum:4}  Weight:  1.00\n")

        f.write("\n//\n\n")

        # Write sequences in blocks of 50 characters
        num_blocks = (nuc_length + 49) // 50
        for block_idx in range(num_blocks):
            start = block_idx * 50
            end = min(start + 50, nuc_length)

            for allele_name, seq in nuc_sequences.items():
                block = seq[start:end]
                # Format in groups of 10
                formatted_block = ' '.join([block[i:i+10] for i in range(0, len(block), 10)])
                f.write(f"{allele_name:<15} {formatted_block}\n")
            f.write("\n")

    print(f"✓ Successfully regenerated {output_path}")
    print(f"  {len(nuc_sequences)} alleles × {nuc_length} bp")

    # Verify against original if it exists
    original_nuc = gen_msf_path.replace('_gen.msf', '_nuc.msf')
    if os.path.exists(original_nuc):
        print(f"\nVerifying regenerated sequences...")
        _, orig_allele_info, orig_sequences = read_msf_file(original_nuc)
        print(f"  Original _nuc.msf: {len(orig_sequences)} alleles, length {len(next(iter(orig_sequences.values())))}")
        print(f"  Regenerated: {len(nuc_sequences)} alleles, length {nuc_length}")

        # Check if biological sequences match for common alleles
        common_alleles = set(nuc_sequences.keys()) & set(orig_sequences.keys())
        if common_alleles:
            # Test a few alleles
            test_count = min(5, len(common_alleles))
            test_alleles = list(common_alleles)[:test_count]
            all_match = True

            for test_allele in test_alleles:
                orig_nogap = orig_sequences[test_allele].replace('.', '').replace('~', '').replace('E', '').replace('N', '')
                regen_nogap = nuc_sequences[test_allele].replace('.', '').replace('~', '').replace('E', '').replace('N', '')
                if orig_nogap != regen_nogap:
                    print(f"  ✗ WARNING: Biological sequences differ for {test_allele}")
                    print(f"    Original: {len(orig_nogap)} bp")
                    print(f"    Regenerated: {len(regen_nogap)} bp")
                    all_match = False
                    break

            if all_match:
                print(f"  ✓ Biological sequences match (tested {test_count} alleles)")

        # Count partial vs full alleles in regenerated
        full_count = len(set(nuc_sequences.keys()) & set(gen_sequences.keys()))
        partial_count = len(nuc_sequences) - full_count
        print(f"  Breakdown: {full_count} full + {partial_count} partial alleles")

    return True


def main():
    # Configuration
    hla_db_dir = "/root/hisat-genotype/indicies/hisatgenotype_db/HLA"
    hla_dat = f"{hla_db_dir}/hla.dat"
    msf_dir = f"{hla_db_dir}/msf"
    output_dir = f"{hla_db_dir}/msf_regenerated"

    # Reference alleles for each gene (from extract_vars output)
    reference_alleles = {
        'A': 'A*03:01:01:01',
        'B': 'B*07:02:01:01',
        'C': 'C*07:02:01:03',
        'DPA1': 'DPA1*01:03:01:02',
        'DPB1': 'DPB1*04:01:01:01',
        'DQB1': 'DQB1*06:02:01:01',
    }

    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    # Process genes (all genes with gap mismatch issues)
    genes_to_process = ['A', 'B', 'C', 'DPA1', 'DPB1', 'DQB1']

    print(f"Regenerating _nuc.msf files...")
    print(f"Input: {msf_dir}")
    print(f"Output: {output_dir}")
    print(f"Reference: {hla_dat}")

    success_count = 0
    for gene in genes_to_process:
        gen_msf = f"{msf_dir}/{gene}_gen.msf"
        output_msf = f"{output_dir}/{gene}_nuc_regenerated.msf"
        ref_allele = reference_alleles.get(gene)

        if not ref_allele:
            print(f"\nSkipping {gene}: no reference allele defined")
            continue

        if not os.path.exists(gen_msf):
            print(f"\nSkipping {gene}: {gen_msf} not found")
            continue

        try:
            if regenerate_nuc_msf(gene, gen_msf, hla_dat, output_msf, ref_allele):
                success_count += 1
        except Exception as e:
            print(f"\nERROR processing {gene}: {e}")
            import traceback
            traceback.print_exc()

    print(f"\n{'='*60}")
    print(f"Regeneration complete: {success_count}/{len(genes_to_process)} genes successful")
    print(f"Output files in: {output_dir}")


if __name__ == "__main__":
    main()
