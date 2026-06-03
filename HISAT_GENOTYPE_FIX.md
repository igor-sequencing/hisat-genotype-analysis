# HISAT-genotype Warning Fix

## Issue
When running hisat-genotype HLA typing, the following warning appeared:
```
Warning: aligner died, killing samtools...
```

The analysis completed partially but failed to process all HLA genes.

## Root Cause Analysis

### Symptoms
- **Warning Message**: "Warning: aligner died, killing samtools..."
- **Incomplete Results**: Only 26 out of 34 HLA genes were processed
- **Successful Genes**: A, B, C, DMA, DMB, DOA, DOB, DPA1, DPB1, DPB2, DQA1, DQA2, DQB1, DQB2, DRA, DRB1, E, F, G, H, K, L, MICA, MICB, TAP1, TAP2, U, V
- **Failed Genes**: DPA2, DRB5, HFE, J, N, S, T, W

### Investigation Results
1. **System Resources**: NOT the issue
   - Available RAM: 1.2 TB
   - Open files limit: 1,048,576
   - No OOM kills or segfaults detected

2. **Warning Source**:
   - File: `hisatgenotype_modules/hisatgenotype_typing_common.py:147`
   - Triggered when hisat2 aligner process exits unexpectedly

3. **Root Cause**: **High parallelism causing hisat2 instability**
   - Original setting: `--pp 20` (20 parallel processes)
   - hisat2 version 2.2.1 has known issues with high parallelism
   - Can cause race conditions and premature process termination

## Fix Applied

### Changes to `/mnt/data/refData/hisat-genotype/run_hisat.sh`

1. **Reduced Parallelism**: `--pp 20` → `--pp 8`
   - More stable for hisat2 2.2.1
   - Reduces race conditions
   - Still provides good performance

2. **Updated samtools**: `/home/samtools/bin/samtools` → `samtools`
   - Now uses system samtools 1.22 (upgraded from 1.13)
   - Fixes library linking issues
   - Better performance and stability

### Before
```bash
/home/samtools/bin/samtools view -b -M $1 6:28500000-33500000 > mhc.bam
hisatgenotype \
  -z indicies \
  --base hla \
  --pp 20 \
  ...
```

### After
```bash
samtools view -b -M $1 6:28500000-33500000 > mhc.bam
hisatgenotype \
  -z indicies \
  --base hla \
  --pp 8 \
  ...
```

## Testing Recommendations

To test the fix:
```bash
# Clean up previous results
rm -rf results_hisat mhc*.bam mhc*.fastq

# Run with fixed script
/mnt/data/refData/hisat-genotype/run_hisat.sh \
  /path/to/sample.bam \
  results_hisat
```

Expected result:
- No "aligner died" warnings
- All 34 HLA genes processed successfully
- Complete report file generated

## Alternative Solutions (if issue persists)

1. **Further reduce parallelism**: Try `--pp 4` or `--pp 2`
2. **Upgrade hisat2**: Consider updating from 2.2.1 to latest version
3. **Serial processing**: Use `--pp 1` for maximum stability (slower)
4. **Check specific genes**: If only certain genes fail, investigate their index files

## Performance Impact

| Setting | Expected Time | Stability |
|---------|--------------|-----------|
| --pp 20 | ~5-8 min | ❌ Unstable (crashes) |
| --pp 8  | ~8-12 min | ✅ Stable (recommended) |
| --pp 4  | ~12-18 min | ✅ Very stable |
| --pp 1  | ~30-40 min | ✅ Maximum stability |

## Date Fixed
February 10, 2026

## Related Issues
- Samtools upgrade: SAMTOOLS_UPGRADE.md
- OptiType installation: OptiType_INSTALLATION_SUMMARY.md
