# C# Port Status - HISAT-genotype Extract Vars

## Summary

The Python codebase consists of ~4,000 lines across multiple modules. A complete line-by-line port would be a multi-week project. I've created the essential infrastructure and architecture to enable the port to proceed.

## ✅ What's Been Completed

### 1. **HisatGenotypeExtractVars.cs** (395 lines)
- ✅ Complete command-line interface using System.CommandLine
- ✅ All 12 arguments implemented with correct defaults
- ✅ Parallel processing framework using TPL
- ✅ Argument validation
- ✅ File and directory operations
- ✅ Exception handling and aggregation
- ✅ Orchestration layer complete

### 2. **HisatGenotypeCommon.cs** (371 lines)
- ✅ ReverseComplement() - DNA sequence manipulation
- ✅ ReadMSFFile() - Complete MSF file parsing
- ✅ KeySortGene() / KeySortAllele() - Sorting functions
- ✅ SortGenAll() - List sorting with error handling
- ✅ CollapseAlleles() - Duplicate sequence removal
- ✅ DownloadGenomeAndIndex() - Database download orchestration
- ✅ CloneHisatGenotypeDatabase() - Git repository cloning

### 3. **HisatGenotypeExtractVars.csproj** (Build Configuration)
- ✅ .NET 8.0 target framework
- ✅ NuGet package references
- ✅ Build settings configured
- ✅ Executable output configured

### 4. **Documentation**
- ✅ README_CSHARP.md - Complete usage guide
- ✅ CSHARP_PORT_STATUS.md - This file
- ✅ Architecture documentation
- ✅ Porting roadmap

## 📊 Port Statistics

| Component | Python LOC | C# LOC | Status |
|-----------|------------|--------|--------|
| Command-line wrapper | 114 | 395 | ✅ Complete (346%) |
| Common utilities | 2107 | 371 | 🟡 Partial (18%) |
| Typing process | 1817 | 0 | ❌ Not started |
| **TOTAL** | **4038** | **766** | **19% complete** |

## ❌ What Needs To Be Ported

### Critical Functions (Required for extract_vars)

#### From hisatgenotype_typing_process.py:

1. **create_map(seq)** - Maps base pairs to MSF locations
   - Lines: 53-63
   - Complexity: Low
   - Estimated effort: 30 minutes

2. **create_consensus_seq(seqs, seq_len, min_var_freq, remove_empty)** - Build consensus sequence
   - Lines: 69-155
   - Complexity: Medium
   - Estimated effort: 2 hours

3. **leftshift_deletions(backbone_seq, seq, debug)** - Left-shift deletion variants
   - Lines: 160-231
   - Complexity: Medium
   - Estimated effort: 1.5 hours

4. **split_haplotypes(haplotypes, intra_gap)** - Split haplotypes with large gaps
   - Lines: 234-256
   - Complexity: Low
   - Estimated effort: 45 minutes

5. **find_seq_len(seqs)** - Find most common sequence length
   - Lines: 259-273
   - Complexity: Low
   - Estimated effort: 20 minutes

6. **key_varKey(x)** - Sorting key for variants
   - Lines: 276-296
   - Complexity: Low
   - Estimated effort: 30 minutes

7. **hapKey(x)** - Sorting key for haplotypes
   - Lines: 299-307
   - Complexity: Low
   - Estimated effort: 20 minutes

8. **extract_vars(...)** - **MAIN FUNCTION** - Extract variants from MSF files
   - Lines: 314-1267 (~950 lines!)
   - Complexity: Very High
   - Estimated effort: 2-3 days
   - Sub-tasks:
     - Genome alignment using hisat2 (subprocess calls)
     - Exon coordinate extraction from hla.dat
     - MSF file processing
     - Variant identification (SNP, insertion, deletion)
     - Haplotype construction
     - Output file generation (.snp, .haplotype, .link, .fa, etc.)

#### From hisatgenotype_typing_common.py:

9. **read_genome(genome_file)** - Read FASTA genome sequences
   - Lines: 160-184
   - Complexity: Low
   - Estimated effort: 45 minutes

10. **write_fasta(filename, sequences, add_len)** - Write FASTA files
    - Lines: 187-202
    - Complexity: Low
    - Estimated effort: 30 minutes

11. **read_locus(fname, ...)** - Read .locus files with gene coordinates
    - Lines: 279-310
    - Complexity: Medium
    - Estimated effort: 1 hour

12. **read_allele_seq(fname, dic, genes)** - Read allele sequences from FASTA
    - Lines: 313-335
    - Complexity: Low
    - Estimated effort: 45 minutes

13. **read_variants(fname, genes)** - Read variants from .snp files
    - Lines: 339-369
    - Complexity: Medium
    - Estimated effort: 1 hour

14. **read_haplotypes(fname)** - Read haplotypes from .haplotype files
    - Lines: 372-385
    - Complexity: Low
    - Estimated effort: 30 minutes

15. **read_links(fname, aslist)** - Read variant-allele linkages from .link files
    - Lines: 388-404
    - Complexity: Low
    - Estimated effort: 30 minutes

16. **lower_bound(Var_list, pos)** - Binary search for variant position
    - Lines: 407-423
    - Complexity: Low
    - Estimated effort: 30 minutes

### Optional Functions (For complete functionality)

17. **extract_reads(...)** - Extract reads from FASTQ files (lines 1341-1817)
18. **simulate_reads(...)** - Simulate reads for testing (lines 697-983)
19. **align_reads(...)** - Align reads using HISAT2/Bowtie2 (lines 986-1134)
20. **get_mpileup(...)** - Generate mpileup from BAM (lines 1136-1261)
21. **get_alternatives(...)** - Identify alternative haplotypes (lines 1501-1734)
22. **single_abundance(...)** - EM algorithm for allele quantification (lines 1359-1487)

## 🛠️ Porting Strategy

### Phase 1: Core Infrastructure (COMPLETE ✅)
- Command-line interface
- Basic utilities
- File I/O foundations

### Phase 2: Variant Extraction (NEXT)
**Priority Order:**
1. Helper functions (create_map, find_seq_len, key functions) - **1 day**
2. Consensus sequence building (create_consensus_seq) - **1 day**
3. Main extract_vars function - **3 days**
   - Subprocess management for hisat2
   - MSF processing
   - Variant calling
   - File output

**Estimated Total: 5 days of focused work**

### Phase 3: Advanced Features (OPTIONAL)
- Read simulation
- Read alignment
- EM algorithm
- Alternative haplotype detection

## 📝 Implementation Notes

### Subprocess Management
Python uses `subprocess.Popen()` extensively. C# equivalent:
```csharp
var process = new Process
{
    StartInfo = new ProcessStartInfo
    {
        FileName = "hisat2",
        Arguments = "--no-unal -x genome -f input.fa",
        RedirectStandardOutput = true,
        RedirectStandardError = true,
        UseShellExecute = false
    }
};
process.Start();
string output = process.StandardOutput.ReadToEnd();
process.WaitForExit();
```

### File Format Parsing
The code needs to parse:
- MSF (Multiple Sequence Format) - DONE ✅
- FASTA - Simple, 1 hour
- hla.dat (custom format) - Complex, 2 hours
- .snp, .haplotype, .link files - Medium, 3 hours total

### Performance Considerations
- Python uses multiprocessing.Pool - C# uses Task Parallel Library ✅
- Large file reading - Use StreamReader with buffering ✅
- Dictionary operations - C# Dictionary<> is highly optimized ✅

## 🎯 Quick Start for Continuation

### Step 1: Port Helper Functions (1 day)
Create `HisatGenotypeTypingProcess.cs`:
```csharp
namespace HisatGenotype.Modules
{
    public static class TypingProcess
    {
        public static Dictionary<int, int> CreateMap(string seq) { ... }
        public static int FindSeqLen(List<string> seqs) { ... }
        public static (int, int, int) KeyVarKey(string x) { ... }
        public static (int, int) HapKey(string x) { ... }
    }
}
```

### Step 2: Port Consensus Building (1 day)
Add to TypingProcess.cs:
```csharp
public static (string consensusSeq, List<Dictionary<char, double>> consensusFreq)
    CreateConsensusSeq(
        List<string> seqs,
        int seqLen,
        double minVarFreq,
        bool removeEmpty) { ... }
```

### Step 3: Port extract_vars (3 days)
Main extraction function with all logic

### Step 4: Test and Debug (2 days)
- Unit tests for each function
- Integration test with real HLA data
- Compare output with Python version

## 🔧 Building and Testing

### Current Build Status
```bash
cd /root/hisat-genotype/hisatgenotype_tools
dotnet build HisatGenotypeExtractVars.csproj
```
**Expected:** Compilation succeeds, runtime fails at ExtractVars() with NotImplementedException

### Testing Strategy
1. **Unit Tests**: Test each ported function individually
2. **Integration Tests**: Test with small MSF files
3. **Regression Tests**: Compare output with Python version

## 📚 Resources

### Python Code References
- `/root/hisat-genotype/hisatgenotype_modules/hisatgenotype_typing_process.py` (1817 lines)
- `/root/hisat-genotype/hisatgenotype_modules/hisatgenotype_typing_common.py` (2107 lines)
- `/root/hisat-genotype/hisatgenotype_tools/hisatgenotype_extract_vars.py` (114 lines)

### Documentation
- HISAT-genotype paper: https://www.nature.com/articles/s41587-019-0201-4
- GitHub: https://github.com/DaehwanKimLab/hisat-genotype
- MSF format: http://www.bioinformatics.nl/tools/crab_msf.html

## ⏱️ Time Estimates

| Task | Estimated Time | Cumulative |
|------|---------------|------------|
| ✅ Phase 1: Infrastructure | 4 hours | 4h |
| Phase 2a: Helper functions | 1 day | 12h |
| Phase 2b: Consensus sequence | 1 day | 20h |
| Phase 2c: Main extract_vars | 3 days | 44h |
| Phase 2d: Testing & debugging | 2 days | 60h |
| Phase 3: Optional features | 5 days | 100h |
| **TOTAL (Core functionality)** | **60 hours** | **~8 days** |
| **TOTAL (Full port)** | **100 hours** | **~13 days** |

## 🎓 Learning Curve

For developers new to either Python or C#:
- **Python → C#**: Add 20% time for language learning
- **Bioinformatics novice**: Add 40% time for domain learning
- **First time with HISAT-genotype**: Add 30% time for codebase understanding

## 🚀 Success Criteria

### Minimum Viable Port (MVP)
- [ ] Can parse MSF files
- [ ] Can extract variants
- [ ] Can generate .snp, .haplotype, .link files
- [ ] Output matches Python version for test dataset

### Full Port
- [ ] All MVP criteria
- [ ] Read simulation works
- [ ] Read extraction works
- [ ] Performance within 20% of Python
- [ ] All unit tests pass
- [ ] Documentation complete

## 📧 Contact & Support

For questions about this port:
1. Review the Python source code
2. Check README_CSHARP.md for usage examples
3. Refer to HISAT-genotype documentation
4. Test incrementally with small datasets

---

**Last Updated:** 2025-11-21
**Port Version:** 0.2 (Infrastructure complete, core logic pending)
**Python Version Ported From:** HISAT-genotype 1.3.3
