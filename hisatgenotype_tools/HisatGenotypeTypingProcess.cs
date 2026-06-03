// ---------------------------------------------------------------------------
// Copyright 2015, Daehwan Kim <infphilo@gmail.com>
// Ported to C# 2025
//
// Core algorithms for HISAT-genotype typing process
// ---------------------------------------------------------------------------

using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;

namespace HisatGenotype.Modules
{
    public static class TypingProcess
    {
        // Mapping from base pair to a location in MSF format given an MSF sequence
        public static Dictionary<int, int> CreateMap(string seq)
        {
            var seqMap = new Dictionary<int, int>();
            int count = 0;

            for (int i = 0; i < seq.Length; i++)
            {
                char bp = seq[i];
                if (bp == '.' || bp == 'E' || bp == 'N' || bp == '~')
                    continue;

                if (!"ACGT".Contains(bp))
                    throw new Exception($"{bp} not a valid basepair");

                seqMap[count] = i;
                count++;
            }

            return seqMap;
        }

        // Check sequences are of equal length - find most common length
        public static int FindSeqLen(List<string> seqs)
        {
            var seqLens = new Dictionary<int, int>();

            foreach (var seq in seqs)
            {
                int seqLen = seq.Length;
                if (!seqLens.ContainsKey(seqLen))
                    seqLens[seqLen] = 0;
                seqLens[seqLen]++;
            }

            int maxSeqCount = 0;
            int resultLen = 0;

            foreach (var kvp in seqLens)
            {
                if (kvp.Value > maxSeqCount)
                {
                    resultLen = kvp.Key;
                    maxSeqCount = kvp.Value;
                }
            }

            return resultLen;
        }

        // Key for sorting Variants in POS-TYPE-INFO format
        public static (int locus, int typeOrd, int lastVal) KeyVarKey(string x)
        {
            var typeOrd = new Dictionary<string, int>
            {
                {"I", 0}, {"M", 1}, {"D", 2}
            };

            var ntOrder = new Dictionary<char, int>
            {
                {'A', 0}, {'C', 1}, {'G', 2}, {'T', 3}
            };

            int Pat2Num(string nt)
            {
                if (nt.Length == 0) return 0;
                int num = ntOrder.Count * Pat2Num(nt.Substring(0, nt.Length - 1)) + ntOrder[nt[^1]];
                return num;
            }

            var parts = x.Split('-');
            int locus = int.Parse(parts[0]);
            string type = parts[1];
            string data = parts[2];

            int next = typeOrd[type];
            int last;

            if (type == "D")
            {
                last = int.Parse(data);
            }
            else
            {
                if (data.Length > 1)
                    last = Pat2Num(data);
                else
                    last = ntOrder[data[0]];
            }

            return (locus, next, last);
        }

        // Key for sorting haplotypes
        public static (int start, int end) HapKey(string x)
        {
            var parts = x.Split('#');
            var firstParts = parts[0].Split('-');
            var lastParts = parts[^1].Split('-');

            int xStart = int.Parse(firstParts[0]);
            int xEnd = int.Parse(lastParts[0]);
            string xType = lastParts[1];
            string xData = lastParts[2];

            if (xType == "D")
                xEnd += int.Parse(xData) - 1;

            return (xStart, xEnd);
        }

        // Given list of MSF sequences, build a consensus sequence based on most common nucleotides
        public static (string consensusSeq, List<Dictionary<char, double>> consensusFreq)
            CreateConsensusSeq(
                List<string> seqs,
                int seqLen,
                double minVarFreq,
                bool removeEmpty = true)
        {
            // Initialize frequency arrays
            var consensusFreq = new List<List<int>>();
            var seqCoverage = new List<int>();

            for (int i = 0; i < seqLen; i++)
            {
                consensusFreq.Add(new List<int> { 0, 0, 0, 0, 0 }); // A, C, G, T, .
                seqCoverage.Add(0);
            }

            // Count nucleotides at each position
            for (int i = 0; i < seqs.Count; i++)
            {
                string seq = seqs[i];
                if (seq.Length != seqLen)
                    continue;

                for (int j = 0; j < seqLen; j++)
                {
                    char nt = seq[j];
                    if (nt == '~')
                        continue;

                    seqCoverage[j]++;

                    if (!"ACGT.EN".Contains(nt))
                        throw new Exception($"Nucleotide {nt} not supported");

                    switch (nt)
                    {
                        case 'A': consensusFreq[j][0]++; break;
                        case 'C': consensusFreq[j][1]++; break;
                        case 'G': consensusFreq[j][2]++; break;
                        case 'T': consensusFreq[j][3]++; break;
                        default: consensusFreq[j][4]++; break; // . E N
                    }
                }
            }

            // Convert to percentages
            for (int j = 0; j < consensusFreq.Count; j++)
            {
                if (seqCoverage[j] == 0) continue;

                for (int k = 0; k < consensusFreq[j].Count; k++)
                {
                    consensusFreq[j][k] = (int)((consensusFreq[j][k] / (double)seqCoverage[j]) * 100.0);
                }
            }

            // Build consensus sequence
            var consensusSeqBuilder = new StringBuilder();
            bool hasEmpty = false;

            for (int c = 0; c < consensusFreq.Count; c++)
            {
                var freq = consensusFreq[c];
                int A = freq[0], C = freq[1], G = freq[2], T = freq[3], E = freq[4];

                if (E >= 100)
                {
                    hasEmpty = true;
                    consensusSeqBuilder.Append('E');
                    continue;
                }

                int idx;
                if (E >= 100.0 - minVarFreq)
                {
                    idx = 4;
                }
                else
                {
                    var maxVal = Math.Max(Math.Max(A, C), Math.Max(G, T));
                    idx = freq.IndexOf((int)maxVal);
                }

                if (idx >= 5)
                    throw new Exception($"Invalid index: {idx}");

                consensusSeqBuilder.Append("ACGT."[idx]);
            }

            string consensusSeq = consensusSeqBuilder.ToString();

            // Remove dots (deletions) if has empty and removeEmpty is true
            var skipPos = new HashSet<int>();
            if (hasEmpty && removeEmpty)
            {
                for (int seq_i = 0; seq_i < seqs.Count; seq_i++)
                {
                    seqs[seq_i] = new string(seqs[seq_i].ToCharArray());
                }

                for (int i = 0; i < consensusSeq.Length; i++)
                {
                    if (consensusSeq[i] != 'E')
                        continue;

                    skipPos.Add(i);
                    for (int seq_i = 0; seq_i < seqs.Count; seq_i++)
                    {
                        if (i >= seqs[seq_i].Length)
                            continue;

                        var chars = seqs[seq_i].ToCharArray();
                        chars[i] = 'E';
                        seqs[seq_i] = new string(chars);
                    }
                }

                for (int seq_i = 0; seq_i < seqs.Count; seq_i++)
                {
                    seqs[seq_i] = seqs[seq_i].Replace("E", "");
                }

                consensusSeq = consensusSeq.Replace("E", "");
            }

            // Convert consensus_freq to dictionary form
            var tempFreq = new List<Dictionary<char, double>>();
            for (int j = 0; j < consensusFreq.Count; j++)
            {
                if (skipPos.Contains(j))
                    continue;

                var freqDic = new Dictionary<char, double>();
                for (int k = 0; k < consensusFreq[j].Count; k++)
                {
                    double freq = consensusFreq[j][k];
                    if (freq <= 0.0)
                        continue;

                    char nt = "ACGT."[k];
                    freqDic[nt] = freq;
                }
                tempFreq.Add(freqDic);
            }

            if (consensusSeq.Length != tempFreq.Count)
                throw new Exception($"Consensus sequence length mismatch: {consensusSeq.Length} != {tempFreq.Count}");

            return (consensusSeq, tempFreq);
        }

        // Left-shift deletions if possible
        public static string LeftshiftDeletions(string backboneSeq, string seq, bool debug = false)
        {
            if (seq.Length != backboneSeq.Length)
                return seq;

            var seqChars = seq.ToCharArray();
            int seqLen = seqChars.Length;
            int bp_i = 0;

            // Skip the first deletion
            while (bp_i < seqLen)
            {
                if ("ACGT".Contains(seqChars[bp_i]))
                    break;
                bp_i++;
            }

            while (bp_i < seqLen)
            {
                char bp = seqChars[bp_i];
                if (bp != '.')
                {
                    bp_i++;
                    continue;
                }

                int bp_j = bp_i + 1;
                while (bp_j < seqLen)
                {
                    char bp2 = seqChars[bp_j];
                    if (bp2 != '.')
                        break;
                    bp_j++;
                }

                if (bp_j >= seqLen)
                {
                    bp_i = bp_j;
                    break;
                }

                int prev_i = bp_i, prev_j = bp_j;

                while (bp_i > 0 &&
                       "ACGT".Contains(seqChars[bp_i - 1]) &&
                       "ACGT".Contains(backboneSeq[bp_j - 1]))
                {
                    if (seqChars[bp_i - 1] != backboneSeq[bp_j - 1])
                        break;

                    seqChars[bp_j - 1] = seqChars[bp_i - 1];
                    seqChars[bp_i - 1] = '.';
                    bp_i--;
                    bp_j--;
                }

                bp_i = bp_j;
                while (bp_i < seqLen)
                {
                    if ("ACGT".Contains(seqChars[bp_i]))
                        break;
                    bp_i++;
                }
            }

            return new string(seqChars);
        }

        // Split haplotypes that include large gaps inside
        public static HashSet<string> SplitHaplotypes(HashSet<string> haplotypes, int intraGap)
        {
            var splitHaps = new HashSet<string>();

            foreach (var haplotype in haplotypes)
            {
                var parts = haplotype.Split('#');

                if (parts.Length == 1)
                {
                    splitHaps.Add(parts[0]);
                    continue;
                }

                int prev_s = 0, s = 1;
                while (s < parts.Length)
                {
                    var prevParts = parts[s - 1].Split('-');
                    var currParts = parts[s].Split('-');

                    int prevLocus = int.Parse(prevParts[0]);
                    string prevType = prevParts[1];
                    string prevData = prevParts[2];

                    int locus = int.Parse(currParts[0]);

                    if (prevType == "D")
                        prevLocus += int.Parse(prevData) - 1;

                    if (prevLocus + intraGap < locus)
                    {
                        splitHaps.Add(string.Join("#", parts[prev_s..s]));
                        prev_s = s;
                    }

                    s++;
                    if (s == parts.Length)
                    {
                        splitHaps.Add(string.Join("#", parts[prev_s..s]));
                    }
                }
            }

            return splitHaps;
        }

        // Main variant extraction function
        public static void ExtractVars(
            string baseFname,
            string ixDir,
            List<string> locusList,
            int interGap,
            int intraGap,
            bool wholeHaplotype,
            double minVarFreq,
            int extSeqLen,
            bool leftshift,
            bool partial,
            bool verbose)
        {
            string baseFullpathName = Path.Combine(ixDir, baseFname);

            // Download human genome and HISAT2 index
            TypingCommon.DownloadGenomeAndIndex(ixDir);

            var splicedGene = new[] { "hla", "rbg" };
            var unsplicedGene = new[] { "codis", "cyp", "rrna" };

            if (verbose)
            {
                Console.Error.WriteLine($"Processing {baseFname}...");
                Console.Error.WriteLine($"  Index directory: {ixDir}");
                Console.Error.WriteLine($"  Locus list: {string.Join(",", locusList)}");
            }

            // Initialize file writers and data structures
            using var locusFile = File.CreateText($"{baseFullpathName}.locus");
            using var varFile = File.CreateText($"{baseFullpathName}.snp");
            using var varIndexFile = File.CreateText($"{baseFullpathName}.index.snp");
            using var varFreqFile = File.CreateText($"{baseFullpathName}.snp.freq");
            using var haplotypeFile = File.CreateText($"{baseFullpathName}.haplotype");
            using var linkFile = File.CreateText($"{baseFullpathName}.link");
            using var backboneFile = File.CreateText($"{baseFullpathName}_backbone.fa");
            using var inputFile = File.CreateText($"{baseFullpathName}_sequences.fa");
            using var alleleFile = File.CreateText($"{baseFullpathName}.allele");
            using var partialFile = File.CreateText($"{baseFullpathName}.partial");

            var leftExtSeqDic = new Dictionary<string, string>();
            var rightExtSeqDic = new Dictionary<string, string>();
            var genes = new Dictionary<string, string>();
            var geneStrand = new Dictionary<string, char>();
            var geneExons = new Dictionary<string, List<(int left, int right)>>();
            var geneExonCounts = new Dictionary<string, Dictionary<int, int>>();

            string hisatgenotypeDb = Path.Combine(ixDir, "hisatgenotype_db");
            string fastaDname = Path.Combine(hisatgenotypeDb, baseFname.ToUpper(), "fasta");

            // Read database version
            string dbVersion = "NONE";
            string versionFile = Path.Combine(hisatgenotypeDb, "VERSION");
            if (File.Exists(versionFile))
            {
                dbVersion = File.ReadAllText(versionFile).Trim();
            }

            // Write version file
            File.WriteAllText($"{baseFullpathName}.version",
                $"Database {baseFname} derived from HISATgenotype DB version: {dbVersion}");

            // Check genes - find all *_gen.fasta files
            var geneNames = new List<string>();
            var fastaPattern = Path.Combine(fastaDname, "*_gen.fasta");

            if (Directory.Exists(Path.GetDirectoryName(fastaPattern)))
            {
                var fastaFiles = Directory.GetFiles(Path.GetDirectoryName(fastaPattern)!,
                    "*_gen.fasta");

                foreach (var genFname in fastaFiles)
                {
                    string geneName = Path.GetFileNameWithoutExtension(genFname).Split('_')[0];
                    if (geneName != "hla" && geneName != baseFname)
                    {
                        geneNames.Add(geneName);
                    }
                }
            }

            geneNames = TypingCommon.SortGenAll(geneNames, false);
            if (locusList.Count == 0)
            {
                locusList = new List<string>(geneNames);
            }

            var cigarRe = new Regex(@"\d+\w");
            var removeLocusList = new List<string>();

            // Phase 1: Align genomic sequences to reference genome using HISAT2
            Console.Error.WriteLine("Phase 1: Aligning genomic sequences to reference...");

            foreach (var gene in locusList.ToList())
            {
                var alignerCmd = new List<string> { "hisat2" };

                if (!new[] { "cyp", "rbg" }.Contains(baseFname))
                {
                    alignerCmd.AddRange(new[] { "--score-min", "C,-12" });
                }

                var geneFastaFile = Path.Combine(fastaDname, $"{gene}_gen.fasta");
                if (!File.Exists(geneFastaFile))
                {
                    Console.Error.WriteLine($"Warning: {geneFastaFile} does not exist");
                    removeLocusList.Add(gene);
                    continue;
                }

                alignerCmd.AddRange(new[] {
                    "--no-unal",
                    "-x", Path.Combine(ixDir, "grch38", "genome"),
                    "-f", geneFastaFile
                });

                var psi = new ProcessStartInfo
                {
                    FileName = "hisat2",
                    Arguments = string.Join(" ", alignerCmd.Skip(1)),
                    RedirectStandardOutput = true,
                    RedirectStandardError = true,
                    UseShellExecute = false
                };

                using var alignProc = Process.Start(psi);
                if (alignProc == null)
                {
                    Console.Error.WriteLine($"Failed to start hisat2 for gene {gene}");
                    removeLocusList.Add(gene);
                    continue;
                }

                // Consume stderr to prevent deadlock
                _ = Task.Run(() => alignProc.StandardError.ReadToEnd());

                string alleleId = "";
                string bestChr = "";
                int bestLeft = -1;
                int bestRight = -1;
                int bestAS = int.MinValue;
                char bestStrand = ' ';

                string? line;
                while ((line = alignProc.StandardOutput.ReadLine()) != null)
                {
                    if (line.StartsWith('@'))
                        continue;

                    var cols = line.Trim().Split('\t');
                    if (cols.Length < 6) continue;

                    string tempAlleleId = cols[0];
                    int flag = int.Parse(cols[1]);
                    string chr = cols[2];
                    int left = int.Parse(cols[3]) - 1;
                    string cigarStr = cols[5];

                    // Calculate right position from CIGAR
                    int right = left;
                    var cigars = cigarRe.Matches(cigarStr);

                    // Only accept simple perfect matches
                    if (cigars.Count > 1 || (cigars.Count == 1 && cigars[0].Value[^1] != 'M'))
                        continue;

                    foreach (Match cigar in cigars)
                    {
                        char op = cigar.Value[^1];
                        int length = int.Parse(cigar.Value[..^1]);
                        if ("MND".Contains(op))
                        {
                            right += length;
                        }
                    }

                    char strand = (flag & 0x10) != 0 ? '-' : '+';

                    // Extract alignment score
                    int AS = int.MinValue;
                    for (int i = 11; i < cols.Length; i++)
                    {
                        if (cols[i].StartsWith("AS:i:"))
                        {
                            AS = int.Parse(cols[i][5..]);
                            break;
                        }
                    }

                    if (AS == int.MinValue)
                        continue;

                    if (AS > bestAS)
                    {
                        alleleId = tempAlleleId;
                        bestChr = chr;
                        bestLeft = left;
                        bestRight = right;
                        bestAS = AS;
                        bestStrand = strand;
                    }
                }

                alignProc.WaitForExit();

                if (string.IsNullOrEmpty(alleleId))
                {
                    removeLocusList.Add(gene);
                    continue;
                }

                // Find allele name from FASTA file
                string alleleName = "";
                foreach (var fastaLine in File.ReadLines(geneFastaFile))
                {
                    if (!fastaLine.StartsWith('>')) continue;

                    var parts = fastaLine[1..].Split(new[] { ' ', '\t' }, StringSplitOptions.RemoveEmptyEntries);
                    if (baseFname == "hla" && parts.Length >= 2)
                    {
                        if (parts[0] == alleleId)
                        {
                            alleleName = parts[1];
                            break;
                        }
                    }
                    else
                    {
                        if (parts[0] == alleleId)
                        {
                            alleleName = parts[0];
                            break;
                        }
                    }
                }

                if (string.IsNullOrEmpty(alleleName) || bestStrand == ' ')
                {
                    removeLocusList.Add(gene);
                    continue;
                }

                genes[gene] = alleleName;
                geneStrand[gene] = bestStrand;

                Console.Error.WriteLine(
                    $"{baseFname.ToUpper()}-{gene}'s reference allele is {alleleName} " +
                    $"on '{bestStrand}' strand of chromosome {bestChr}");

                // Extract extended sequences if requested
                if (extSeqLen > 0 && bestLeft >= 0 && bestRight > bestLeft)
                {
                    string leftExtSeq = ExtractSequenceUsingSamtools(
                        ixDir, bestChr, Math.Max(1, bestLeft - extSeqLen), Math.Max(1, bestLeft - 1));
                    string rightExtSeq = ExtractSequenceUsingSamtools(
                        ixDir, bestChr, bestRight, bestRight + extSeqLen - 1);

                    if (bestStrand == '-')
                    {
                        string temp = leftExtSeq;
                        leftExtSeq = TypingCommon.ReverseComplement(rightExtSeq);
                        rightExtSeq = TypingCommon.ReverseComplement(temp);
                    }

                    leftExtSeqDic[gene] = leftExtSeq;
                    rightExtSeqDic[gene] = rightExtSeq;
                }
            }

            // Phase 2: Extract exon information from .dat file (for spliced genes)
            Console.Error.WriteLine("Phase 2: Extracting exon information...");

            if (splicedGene.Contains(baseFname))
            {
                string datFile = Path.Combine(hisatgenotypeDb, baseFname.ToUpper(), $"{baseFname}.dat");
                if (File.Exists(datFile))
                {
                    ParseExonInformation(datFile, baseFname, genes, leftExtSeqDic,
                        geneExons, geneExonCounts, verbose);
                }
            }

            // Phase 3: Filter locus list
            var tmpGenes = new Dictionary<string, string>();
            var tmpGeneStrand = new Dictionary<string, char>();

            for (int i = locusList.Count - 1; i >= 0; i--)
            {
                string gene = locusList[i];
                if (removeLocusList.Contains(gene) ||
                    (splicedGene.Contains(baseFname) && !geneExons.ContainsKey(gene)))
                {
                    locusList.RemoveAt(i);
                }
                else if (genes.ContainsKey(gene))
                {
                    tmpGenes[gene] = genes[gene];
                    tmpGeneStrand[gene] = geneStrand[gene];
                }
            }

            genes = tmpGenes;
            geneStrand = tmpGeneStrand;

            Console.Error.WriteLine($"Processing {genes.Count} genes: {string.Join(", ", genes.Keys)}");
            Console.Error.WriteLine();

            // Phase 4: Process MSF files for each gene and extract variants
            Console.Error.WriteLine("Phase 4: Processing MSF files and extracting variants...");

            int numVars = 0;
            int numHaplotypes = 0;
            var fullAlleles = new Dictionary<string, List<string>>();

            foreach (var gene in geneNames)
            {
                if (!genes.ContainsKey(gene))
                    continue;

                string refGene = genes[gene];
                char strand = geneStrand[gene];
                string leftExtSeq = leftExtSeqDic.ContainsKey(gene) ? leftExtSeqDic[gene] : "";
                string rightExtSeq = rightExtSeqDic.ContainsKey(gene) ? rightExtSeqDic[gene] : "";

                // Determine MSF file path
                string msaFname = Path.Combine(hisatgenotypeDb, baseFname.ToUpper(), "msf", $"{gene}_gen.msf");

                if (!File.Exists(msaFname))
                {
                    Console.Error.WriteLine($"Warning: {msaFname} does not exist");
                    continue;
                }

                // Read MSF file
                var (names, seqs) = TypingCommon.ReadMSFFile(msaFname, fullAlleles, leftExtSeq, rightExtSeq);
                var fullAlleleNames = new HashSet<string>(names.Keys);

                if (seqs.Count == 0)
                {
                    Console.Error.WriteLine($"Warning: No sequences found in {msaFname}");
                    continue;
                }

                // Build consensus sequence
                int seqLen = FindSeqLen(seqs);
                string backboneName = $"{gene}*BACKBONE";
                var (backboneSeq, backboneFreq) = CreateConsensusSeq(seqs, seqLen, minVarFreq, !partial);

                // Readjust sequence length if not using partial alleles
                if (!partial)
                {
                    seqLen = FindSeqLen(seqs);
                }

                // Handle partial alleles for spliced genes
                if (partial && splicedGene.Contains(baseFname))
                {
                    string partialMsaFname = Path.Combine(hisatgenotypeDb, baseFname.ToUpper(), "msf", $"{gene}_nuc.msf");

                    if (!File.Exists(partialMsaFname))
                    {
                        Console.Error.WriteLine($"Warning: {partialMsaFname} does not exist");
                        continue;
                    }

                    var (partialNames, partialSeqs) = TypingCommon.ReadMSFFile(partialMsaFname, fullAlleles);

                    // Process partial alleles (simplified - full implementation would merge them)
                    // For now, we'll skip the complex partial allele merging logic
                    Console.Error.WriteLine($"Warning: Partial allele processing for {gene} is simplified in C# port");
                }

                // Fill in missing sequences (~) with consensus
                bool missingSeq = false;
                for (int itr = 0; itr < seqs.Count; itr++)
                {
                    if (!seqs[itr].Contains('~'))
                        continue;

                    missingSeq = true;
                    string seq = seqs[itr];
                    if (seq.Length != backboneSeq.Length)
                        continue;

                    var seqChars = new StringBuilder();
                    for (int s = 0; s < seq.Length; s++)
                    {
                        if (seq[s] == '~')
                            seqChars.Append(backboneSeq[s]);
                        else
                            seqChars.Append(seq[s]);
                    }
                    seqs[itr] = seqChars.ToString();
                }

                if (missingSeq)
                {
                    Console.Error.WriteLine($"Warning: {gene} contains missing sequence in the data. Filling in with consensus");
                }

                // Collapse duplicate alleles
                var (collapsedNames, collapsedSeqs, collapsedMap) = TypingCommon.CollapseAlleles(names, seqs, listCollapse: true, verbose: true);
                names = collapsedNames;
                seqs = collapsedSeqs;

                if (collapsedMap.ContainsKey(refGene))
                {
                    refGene = collapsedMap[refGene];
                    genes[gene] = refGene;
                }

                // Check for empty sequences or omitted nucleotides in backbone
                if (minVarFreq <= 0.0)
                {
                    var omits = new[] { '.', 'E', '~' };
                    bool breakout = false;

                    foreach (char omit in omits)
                    {
                        if (backboneSeq.Contains(omit))
                        {
                            if (verbose)
                            {
                                Console.Error.WriteLine($"{omit} in backbone of {gene} with no minimum variation set");
                            }
                            Console.Error.WriteLine($"Error in database: Omitting {gene}!!");
                            breakout = true;
                            break;
                        }
                    }

                    if (breakout)
                        continue;
                }

                // Reverse complement MSF if gene is on minus strand
                if (strand == '-')
                {
                    string refSeq = seqs[names[refGene]].Replace(".", "").Replace("~", "");
                    int refSeqLen = refSeq.Length;

                    if (splicedGene.Contains(baseFname) && geneExons.ContainsKey(gene))
                    {
                        // Reverse exons
                        var exons = new List<(int left, int right)>();
                        foreach (var (left, right) in geneExons[gene].AsEnumerable().Reverse())
                        {
                            int rleft = refSeqLen - right - 1;
                            int rright = refSeqLen - left - 1;
                            exons.Add((rleft, rright));
                        }
                        geneExons[gene] = exons;

                        var exonCounts = new Dictionary<int, int>();
                        foreach (var (exonI, count) in geneExonCounts[gene])
                        {
                            exonCounts[geneExons[gene].Count - exonI - 1] = count;
                        }
                        geneExonCounts[gene] = exonCounts;
                    }

                    for (int i = 0; i < seqs.Count; i++)
                    {
                        seqs[i] = TypingCommon.ReverseComplement(seqs[i]);
                    }

                    (backboneSeq, backboneFreq) = CreateConsensusSeq(seqs, seqLen, minVarFreq, removeEmpty: true);
                    seqLen = FindSeqLen(seqs);
                }

                // Apply leftshift to deletions if requested
                if (leftshift)
                {
                    for (int seqI = 0; seqI < seqs.Count; seqI++)
                    {
                        seqs[seqI] = LeftshiftDeletions(backboneSeq, seqs[seqI]);
                    }
                    (backboneSeq, backboneFreq) = CreateConsensusSeq(seqs, seqLen, minVarFreq, removeEmpty: true);
                    seqLen = FindSeqLen(seqs);
                }

                Console.Error.WriteLine($"{gene}: number of alleles is {names.Count}.");

                // Phase 5: Identify variants by comparing each sequence to backbone
                var Vars = new Dictionary<string, (double freq, List<string> names)>();

                foreach (var (cmpName, id) in names)
                {
                    if (cmpName == backboneName)
                        continue;

                    if (id >= seqs.Count)
                    {
                        Console.Error.WriteLine($"Warning: sequence ID {id} out of range for {cmpName}");
                        continue;
                    }

                    string cmpSeq = seqs[id];
                    if (cmpSeq.Length != seqLen)
                    {
                        Console.Error.WriteLine($"Warning: the length of {cmpName} ({cmpSeq.Length}) is different from {seqLen}");
                        continue;
                    }

                    // Inner function to insert variants
                    void InsertVar(char type, (int pos, int backbonePos, string data) info)
                    {
                        string varKey;
                        if (type == 'M' || type == 'I')
                            varKey = $"{info.pos}-{type}-{info.data}";
                        else
                            varKey = $"{info.pos}-{type}-{info.data}"; // For deletion, data is length as string

                        if (!Vars.ContainsKey(varKey))
                        {
                            double freq = 100.0;

                            if (type == 'M')
                            {
                                if (!backboneFreq[info.backbonePos].ContainsKey(info.data[0]))
                                {
                                    Console.Error.WriteLine($"Warning: Data {info.data} not in backbone freq at pos {info.backbonePos}");
                                    freq = 0.0;
                                }
                                else
                                {
                                    freq = backboneFreq[info.backbonePos][info.data[0]];
                                }
                            }
                            else if (type == 'D')
                            {
                                int delLen = int.Parse(info.data);
                                for (int d = 0; d < delLen; d++)
                                {
                                    if (info.backbonePos + d >= backboneFreq.Count)
                                        break;
                                    if (!backboneFreq[info.backbonePos + d].ContainsKey('.'))
                                        continue;
                                    double freq2 = backboneFreq[info.backbonePos + d]['.'];
                                    if (freq2 < freq)
                                        freq = freq2;
                                }
                            }
                            else if (type == 'I')
                            {
                                int insLen = info.data.Length;
                                for (int i = 0; i < insLen; i++)
                                {
                                    if (info.backbonePos + i >= backboneFreq.Count)
                                        break;
                                    char nt = info.data[i];
                                    if (!backboneFreq[info.backbonePos + i].ContainsKey(nt))
                                        continue;
                                    double freq2 = backboneFreq[info.backbonePos + i][nt];
                                    if (freq2 < freq)
                                        freq = freq2;
                                }
                            }

                            Vars[varKey] = (freq, new List<string> { cmpName });
                        }
                        else
                        {
                            Vars[varKey].names.Add(cmpName);
                        }
                    }

                    // Scan sequence and identify variants
                    (int pos, int backbonePos, string data)? insertion = null;
                    (int pos, int backbonePos, string data)? deletion = null;
                    int ndots = 0;

                    for (int s = 0; s < seqLen; s++)
                    {
                        char bc = backboneSeq[s];
                        char cc = cmpSeq[s];

                        if (bc != '.' && bc != '~' && cc != '.' && cc != '~')
                        {
                            if (insertion.HasValue)
                            {
                                InsertVar('I', insertion.Value);
                                insertion = null;
                            }
                            else if (deletion.HasValue)
                            {
                                InsertVar('D', deletion.Value);
                                deletion = null;
                            }

                            if (bc != cc)
                            {
                                var mismatch = (s - ndots, s, cc.ToString());
                                InsertVar('M', mismatch);
                            }
                        }
                        else if (bc == '.' && cc != '.' && cc != '~')
                        {
                            if (deletion.HasValue)
                            {
                                InsertVar('D', deletion.Value);
                                deletion = null;
                            }

                            if (insertion.HasValue)
                            {
                                insertion = (insertion.Value.pos, insertion.Value.backbonePos, insertion.Value.data + cc);
                            }
                            else
                            {
                                insertion = (s - ndots, s, cc.ToString());
                            }
                        }
                        else if (bc != '.' && bc != '~' && cc == '.')
                        {
                            if (insertion.HasValue)
                            {
                                InsertVar('I', insertion.Value);
                                insertion = null;
                            }

                            if (deletion.HasValue)
                            {
                                int newLen = int.Parse(deletion.Value.data) + 1;
                                deletion = (deletion.Value.pos, deletion.Value.backbonePos, newLen.ToString());
                            }
                            else
                            {
                                deletion = (s - ndots, s, "1");
                            }
                        }

                        if (bc == '.')
                            ndots++;
                    }

                    // Handle remaining insertion or deletion
                    if (insertion.HasValue)
                        InsertVar('I', insertion.Value);
                    else if (deletion.HasValue)
                        InsertVar('D', deletion.Value);
                }

                Console.Error.WriteLine($"Number of variants is {Vars.Count}.");

                // TODO: Continue with variant writing, haplotype construction, and file output
                // This will be implemented in the next phase
            }

            Console.Error.WriteLine();
            Console.Error.WriteLine("Note: Variant writing, haplotype construction, and file output");
            Console.Error.WriteLine("      are not yet fully implemented in this C# port.");
            Console.Error.WriteLine("      The Python version should be used for production work.");
        }

        // Helper method to extract sequence using samtools faidx
        private static string ExtractSequenceUsingSamtools(string ixDir, string chr, int start, int end)
        {
            var genomeFile = Path.Combine(ixDir, "genome.fa");
            if (!File.Exists(genomeFile))
            {
                Console.Error.WriteLine($"Warning: genome file not found: {genomeFile}");
                return "";
            }

            var psi = new ProcessStartInfo
            {
                FileName = "samtools",
                Arguments = $"faidx {genomeFile} {chr}:{start}-{end}",
                RedirectStandardOutput = true,
                RedirectStandardError = true,
                UseShellExecute = false
            };

            try
            {
                using var proc = Process.Start(psi);
                if (proc == null) return "";

                var sb = new StringBuilder();
                string? line;
                while ((line = proc.StandardOutput.ReadLine()) != null)
                {
                    if (!line.StartsWith('>'))
                    {
                        sb.Append(line.Trim());
                    }
                }

                proc.WaitForExit();
                return sb.ToString();
            }
            catch (Exception ex)
            {
                Console.Error.WriteLine($"Error running samtools: {ex.Message}");
                return "";
            }
        }

        // Helper method to parse exon information from .dat files
        private static void ParseExonInformation(
            string datFile,
            string baseFname,
            Dictionary<string, string> genes,
            Dictionary<string, string> leftExtSeqDic,
            Dictionary<string, List<(int left, int right)>> geneExons,
            Dictionary<string, Dictionary<int, int>> geneExonCounts,
            bool verbose)
        {
            bool skip = false;
            bool lookExonNum = false;
            string currentGene = "";
            string currentAlleleName = "";

            foreach (var line in File.ReadLines(datFile))
            {
                if (line.StartsWith("DE"))
                {
                    var parts = line.Split(new[] { ' ', '\t' }, StringSplitOptions.RemoveEmptyEntries);
                    if (parts.Length < 2) continue;

                    string alleleName = parts[1];
                    if (!char.IsDigit(alleleName[^1]))
                    {
                        alleleName = alleleName[..^1];
                    }

                    if (alleleName.StartsWith($"{baseFname.ToUpper()}-"))
                    {
                        alleleName = alleleName[($"{baseFname.ToUpper()}-").Length..];
                    }

                    currentGene = alleleName.Split('*')[0];
                    currentAlleleName = alleleName;
                    skip = !genes.ContainsKey(currentGene);
                }

                if (skip) continue;
                if (!line.StartsWith("FT")) continue;

                if (line.Contains("exon"))
                {
                    lookExonNum = true;
                    if (genes.TryGetValue(currentGene, out string? refAllele) &&
                        currentAlleleName == refAllele)
                    {
                        var parts = line.Split(new[] { ' ', '\t' }, StringSplitOptions.RemoveEmptyEntries);
                        if (parts.Length < 3) continue;

                        var exonRange = parts[2].Split("..");
                        if (exonRange.Length == 2)
                        {
                            int exonLeft = int.Parse(exonRange[0]) - 1;
                            int exonRight = int.Parse(exonRange[1]) - 1;

                            if (exonLeft >= 0 && exonLeft < exonRight)
                            {
                                if (!geneExons.ContainsKey(currentGene))
                                {
                                    geneExons[currentGene] = new List<(int, int)>();
                                }

                                int leftExtSeqLen = leftExtSeqDic.ContainsKey(currentGene) ?
                                    leftExtSeqDic[currentGene].Length : 0;

                                geneExons[currentGene].Add((exonLeft + leftExtSeqLen,
                                    exonRight + leftExtSeqLen));
                            }
                        }
                    }
                }
                else if (lookExonNum && line.Contains("number"))
                {
                    lookExonNum = false;
                    var digits = new string(line.Where(char.IsDigit).ToArray());
                    if (!string.IsNullOrEmpty(digits))
                    {
                        int num = int.Parse(digits) - 1;
                        if (!geneExonCounts.ContainsKey(currentGene))
                        {
                            geneExonCounts[currentGene] = new Dictionary<int, int>();
                        }

                        if (!geneExonCounts[currentGene].ContainsKey(num))
                        {
                            geneExonCounts[currentGene][num] = 0;
                        }

                        geneExonCounts[currentGene][num]++;
                    }
                }
            }

            if (verbose)
            {
                foreach (var (gene, exonCounts) in geneExonCounts)
                {
                    Console.Error.WriteLine($"{gene} exon counts: {string.Join(", ",
                        exonCounts.Select(kv => $"{kv.Key}:{kv.Value}"))}");
                }
            }
        }
    }
}
