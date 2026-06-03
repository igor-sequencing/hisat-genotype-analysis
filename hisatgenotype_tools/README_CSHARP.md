# C# Port of hisatgenotype_extract_vars.py

This directory contains a C# port of the `hisatgenotype_extract_vars.py` script.

## Files

- **HisatGenotypeExtractVars.cs** - Main C# implementation
- **HisatGenotypeExtractVars.csproj** - .NET project file
- **README_CSHARP.md** - This file

## Overview

The C# version maintains the same command-line interface and functionality as the original Python script:

### Features Implemented

✅ **Command-line argument parsing** using System.CommandLine
- All original options supported (--base, --locus-list, -z, --inter-gap, --intra-gap, etc.)
- Same default values as Python version
- Validation of arguments (e.g., inter-gap < intra-gap)

✅ **Parallel processing** using Task Parallel Library (TPL)
- Replaces Python's multiprocessing.Pool
- Configurable thread count (-p option)
- Exception handling and aggregation

✅ **File and directory operations**
- Index directory resolution
- Database path checking
- Support for hg_ix.link file

✅ **Git integration** placeholder
- Clone hisatgenotype_db repository when needed

### Features NOT Implemented (Require Additional Porting)

❌ **hisatgenotype_typing_process module**
- The core `extract_vars` function needs to be ported from Python
- File: `hisatgenotype_modules/hisatgenotype_typing_process.py`
- This is ~1000+ lines of MSF parsing, variant extraction logic

❌ **hisatgenotype_typing_common module**
- Various utility functions
- File: `hisatgenotype_modules/hisatgenotype_typing_common.py`

## Requirements

- .NET 8.0 SDK or later
- System.CommandLine NuGet package (automatically restored)

## Building

```bash
cd /root/hisat-genotype/hisatgenotype_tools
dotnet build HisatGenotypeExtractVars.csproj
```

## Running

```bash
# After building:
dotnet run --project HisatGenotypeExtractVars.csproj -- --base hla -p 4 -v -z ../indicies

# Or after publishing:
./bin/Release/net8.0/hisatgenotype_extract_vars --base hla -p 4 -v -z ../indicies
```

## Usage

The C# version supports all the same options as the Python version:

```bash
hisatgenotype_extract_vars [OPTIONS]

Options:
  --base, --base-fname <base>      Base file name (e.g., hla, rbg, codis)
  --locus-list <list>              Comma-separated gene names
  -z, --index-dir <dir>            Index directory location
  --inter-gap <n>                  Max distance for variants in same haplotype (default: 30)
  --intra-gap <n>                  Break haplotype threshold (default: 50)
  --whole-haplotype                Include whole haplotypes
  --min-var-freq <freq>            Minimum variant frequency (default: 0.0)
  --ext-seq <len>                  Extended sequence length (default: 0)
  --leftshift                      Shift deletions to leftmost
  --no-partial                     Exclude partial alleles
  -p, --threads <n>                Number of threads (default: 1)
  -v, --verbose                    Print statistics to stderr
  --help                           Show help
  --version                        Show version
```

## Architecture Differences

### Python vs C# Implementation

| Aspect | Python | C# |
|--------|--------|-----|
| **Argument Parsing** | argparse | System.CommandLine |
| **Parallel Processing** | multiprocessing.Pool | Task Parallel Library (TPL) |
| **Synchronization** | multiprocessing.Lock | SemaphoreSlim |
| **File I/O** | os, sys | System.IO |
| **Process Execution** | subprocess | System.Diagnostics.Process |
| **Exception Handling** | try/except | try/catch with AggregateException |

## What Needs to Be Done for Full Functionality

To make this a complete, functional C# port, you need to:

1. **Port hisatgenotype_typing_process.py** (~1000+ lines)
   - MSF file parsing
   - Variant extraction algorithms
   - Haplotype processing
   - FASTA file generation

2. **Port hisatgenotype_typing_common.py**
   - Sequence alignment utilities
   - Database cloning/management
   - Common data structures

3. **Port hisatgenotype_args.py** (partially done)
   - All argument definitions are implemented
   - Some validation logic may need porting

4. **Testing**
   - Unit tests for all components
   - Integration tests with real HLA data
   - Performance benchmarking vs Python version

## Performance Considerations

**Expected Benefits of C# Port:**
- ✅ Faster startup time (no Python interpreter)
- ✅ Better memory management (compiled vs interpreted)
- ✅ Native threading (vs GIL limitations in Python)
- ✅ Potential for SIMD optimizations
- ✅ AOT compilation option (Native AOT)

**Trade-offs:**
- ⚠️ Larger binary size
- ⚠️ More complex deployment (need .NET runtime)
- ⚠️ String manipulation might be slower than Python for some operations

## Development Notes

### Current Status
- ✅ Wrapper/orchestration layer complete
- ✅ Command-line interface complete
- ✅ Parallel processing framework complete
- ❌ Core variant extraction logic needs porting
- ❌ MSF file parsing needs porting
- ❌ Database utilities need porting

### Code Organization
The C# version follows .NET conventions:
- PascalCase for public members
- camelCase for parameters and local variables
- Async/await for I/O operations
- LINQ for collection operations
- Nullable reference types enabled

### Testing Against Python Version

To verify compatibility:

```bash
# Run Python version
python hisatgenotype_extract_vars.py --base hla -p 4 -v -z ../indicies

# Run C# version (once core logic is ported)
dotnet run -- --base hla -p 4 -v -z ../indicies

# Compare outputs
diff python_output/ csharp_output/
```

## Contributing

If you continue this port, focus on:

1. **hisatgenotype_typing_process.ExtractVars()** - Core function
2. **MSF file parsing** - Read/write multiple sequence alignments
3. **Variant calling logic** - SNP, insertion, deletion detection
4. **FASTA generation** - Output file creation

## License

Same as original: GNU General Public License v3.0 or later

## Original Authors

- Daehwan Kim <infphilo@gmail.com>
- C# port: 2025

## References

- Original Python code: `/root/hisat-genotype/hisatgenotype_tools/hisatgenotype_extract_vars.py`
- HISAT-genotype repository: https://github.com/DaehwanKimLab/hisat-genotype
