#!/usr/bin/env python3
import sys
sys.path.insert(0, 'hisatgenotype_modules')

from hisatgenotype_typing_process import extract_vars

# Call extract_vars directly without multiprocessing
extract_vars(
    base_fname='hla',
    ix_dir='indicies',
    locus_list=[],  # Empty list means all genes
    inter_gap=30,
    intra_gap=50,
    whole_haplotype=True,  # Changed to True for comparison
    min_var_freq=0.0,
    ext_seq_len=0,
    leftshift=False,
    partial=True,
    verbose=True
)

print("Extract vars completed successfully!")
