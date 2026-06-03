#!/usr/bin/env python3
"""
Extract HLA region (chr6:28477797-33448354) from reference genome
"""

import sys

def read_fasta_index(fai_file):
    """Read FASTA index file to get chromosome positions"""
    index = {}
    with open(fai_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            chrom = fields[0]
            length = int(fields[1])
            offset = int(fields[2])
            line_bases = int(fields[3])
            line_width = int(fields[4])
            index[chrom] = {
                'length': length,
                'offset': offset,
                'line_bases': line_bases,
                'line_width': line_width
            }
    return index

def extract_region(fasta_file, chrom, start, end, index):
    """Extract specific region from FASTA file"""
    if chrom not in index:
        raise ValueError(f"Chromosome {chrom} not found in index")

    info = index[chrom]

    # Convert to 0-based coordinates
    start_0 = start - 1
    end_0 = end

    # Calculate file position for start
    lines_before_start = start_0 // info['line_bases']
    pos_in_line = start_0 % info['line_bases']
    file_start = info['offset'] + lines_before_start * info['line_width'] + pos_in_line

    # Read the sequence
    sequence = []
    with open(fasta_file, 'rb') as f:
        f.seek(file_start)
        bases_to_read = end_0 - start_0
        bases_read = 0

        while bases_read < bases_to_read:
            chunk = f.read(8192).decode('ascii')
            if not chunk:
                break
            # Remove newlines
            chunk_clean = chunk.replace('\n', '').replace('\r', '')
            sequence.append(chunk_clean)
            bases_read += len(chunk_clean)

    full_sequence = ''.join(sequence)[:bases_to_read]
    return full_sequence

def write_fasta(output_file, header, sequence, line_width=70):
    """Write sequence to FASTA file"""
    with open(output_file, 'w') as f:
        f.write(f">{header}\n")
        for i in range(0, len(sequence), line_width):
            f.write(sequence[i:i+line_width] + '\n')

def main():
    # Parameters
    ref_fasta = '/refData/1065/1065.fa'
    ref_index = ref_fasta + '.fai'
    chrom = '6'
    start = 28477797
    end = 33448354
    output_file = 'HLA.fasta'

    print(f"Reading FASTA index from {ref_index}...")
    index = read_fasta_index(ref_index)

    print(f"Extracting {chrom}:{start}-{end}...")
    sequence = extract_region(ref_fasta, chrom, start, end, index)

    print(f"Extracted {len(sequence):,} bases")

    header = f"{chrom}:{start}-{end}"
    print(f"Writing to {output_file}...")
    write_fasta(output_file, header, sequence)

    print(f"Done! HLA region saved to {output_file}")
    print(f"File size: {len(sequence):,} bases")

if __name__ == '__main__':
    main()
