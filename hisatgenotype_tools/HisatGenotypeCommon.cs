// ---------------------------------------------------------------------------
// Copyright 2017, Daehwan Kim <infphilo@gmail.com>
// Ported to C# 2025
//
// Common utility functions for HISAT-genotype
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
    public static class TypingCommon
    {
        // Sequence processing
        public static string ReverseComplement(string seq)
        {
            var compTable = new Dictionary<char, char>
            {
                {'A', 'T'}, {'C', 'G'}, {'G', 'C'}, {'T', 'A'}
            };

            var sb = new StringBuilder(seq.Length);
            for (int i = seq.Length - 1; i >= 0; i--)
            {
                char c = seq[i];
                sb.Append(compTable.ContainsKey(c) ? compTable[c] : c);
            }
            return sb.ToString();
        }

        // MSF file reading
        public static (Dictionary<string, int> names, List<string> seqs) ReadMSFFile(
            string fname,
            Dictionary<string, List<string>> fullAlleles,
            string leftExtSeq = "",
            string rightExtSeq = "")
        {
            var names = new Dictionary<string, int>();
            var seqs = new List<string>();

            using (var reader = new StreamReader(fname))
            {
                string line;
                while ((line = reader.ReadLine()) != null)
                {
                    line = line.Trim();

                    if (string.IsNullOrWhiteSpace(line) ||
                        !char.IsLetterOrDigit(line[0]) ||
                        line.StartsWith("MSF") ||
                        line.StartsWith("PileUp"))
                    {
                        continue;
                    }

                    if (line.StartsWith("Name"))
                    {
                        try
                        {
                            var parts = line.Split('\t', StringSplitOptions.RemoveEmptyEntries);
                            if (parts.Length == 0) continue;

                            var nameParts = parts[0].Split(new[] {' '}, StringSplitOptions.RemoveEmptyEntries);
                            if (nameParts.Length < 2) continue;

                            string name = nameParts[1];

                            if (names.ContainsKey(name))
                            {
                                Console.Error.WriteLine($"Warning: {name} is found more than once in Names");
                                continue;
                            }

                            names[name] = names.Count;
                        }
                        catch
                        {
                            continue;
                        }
                    }
                    else
                    {
                        if (seqs.Count == 0)
                        {
                            for (int i = 0; i < names.Count; i++)
                            {
                                seqs.Add(leftExtSeq);
                            }
                        }

                        try
                        {
                            var cols = line.Split(new[] {' ', '\t'}, StringSplitOptions.RemoveEmptyEntries);
                            if (cols.Length == 0) continue;

                            string name = cols[0];
                            var fives = cols.Skip(1).ToList();

                            if (fives.Count == 0) continue;

                            if (!names.ContainsKey(name))
                            {
                                names[name] = names.Count;
                            }

                            int id = names[name];
                            if (id >= seqs.Count)
                            {
                                if (id != seqs.Count) throw new Exception($"Unexpected sequence ID: {id}");
                                seqs.Add(leftExtSeq);
                            }

                            seqs[id] += string.Join("", fives);

                            // Add sub-names of the allele
                            var groups = name.Split(':');
                            string subName = "";
                            for (int i = 0; i < groups.Length - 1; i++)
                            {
                                if (subName != "") subName += ":";
                                subName += groups[i];

                                if (!fullAlleles.ContainsKey(subName))
                                {
                                    fullAlleles[subName] = new List<string> { name };
                                }
                                else
                                {
                                    fullAlleles[subName].Add(name);
                                }
                            }
                        }
                        catch
                        {
                            continue;
                        }
                    }
                }
            }

            if (rightExtSeq.Length > 0)
            {
                for (int i = 0; i < seqs.Count; i++)
                {
                    seqs[i] += rightExtSeq;
                }
            }

            return (names, seqs);
        }

        // Sorting functions
        public static (string chars, int nums) KeySortGene(string x)
        {
            var digits = new List<char>();
            var chars = new List<char>();

            foreach (char y in x)
            {
                if (char.IsDigit(y))
                    digits.Add(y);
                else
                    chars.Add(y);
            }

            if (digits.Count == 0) digits.Add('-'); digits.Add('1');
            if (chars.Count == 0) chars.Add('\0');

            int nums = int.Parse(new string(digits.ToArray()));
            string strs = new string(chars.ToArray());

            return (strs, nums);
        }

        public static List<string> SortGenAll(List<string> list_, bool alleles = false)
        {
            try
            {
                if (alleles)
                {
                    return list_.OrderBy(x => KeySortAllele(x)).ToList();
                }
                else
                {
                    return list_.OrderBy(x => KeySortGene(x)).ToList();
                }
            }
            catch (Exception ex)
            {
                Console.Error.WriteLine("Error in sorting list of alleles or genes!!!");
                Console.Error.WriteLine(ex.Message);
                Environment.Exit(1);
                return null;
            }
        }

        private static (string, int, int, int, int, int) KeySortAllele(string x)
        {
            var parts = x.Split('*');
            string gene = parts[0];
            string allele = parts[1];

            var (gen, val) = KeySortGene(gene);

            var alleleFields = allele.Split(':')
                .Select(f => int.Parse(Regex.Replace(f, @"[^0-9]", "")))
                .ToList();

            while (alleleFields.Count < 4)
            {
                alleleFields.Add(-1);
            }

            return (gen, val, alleleFields[0], alleleFields[1], alleleFields[2], alleleFields[3]);
        }

        // Collapse duplicate alleles
        public static (Dictionary<string, int> index, List<string> seqs, Dictionary<string, string> collapsed)
            CollapseAlleles(
                Dictionary<string, int> index,
                List<string> seqs,
                string emptySeq = "",
                bool listCollapse = false,
                bool verbose = false)
        {
            var remove = new List<(int, string)>();
            var colIndex = new Dictionary<string, string>();

            foreach (var kvp_i in index)
            {
                string allele_i = kvp_i.Key;
                int index_i = kvp_i.Value;

                if (remove.Any(x => x.Item2 == allele_i)) continue;

                string seq_i = seqs[index_i];
                string seq_i_strip = seq_i.Replace(emptySeq, "").Replace(".", "");

                foreach (var kvp_j in index)
                {
                    string allele_j = kvp_j.Key;
                    int index_j = kvp_j.Value;

                    if (remove.Any(x => x.Item2 == allele_j)) continue;

                    string seq_j = seqs[index_j];
                    string seq_j_strip = seq_j.Replace(emptySeq, "").Replace(".", "");

                    if (allele_i == allele_j) continue;

                    if (seq_i == seq_j && string.Compare(allele_i, allele_j) > 0)
                    {
                        if (allele_i.Length <= allele_j.Length)
                        {
                            if (verbose)
                                Console.Error.WriteLine($"\t\t {allele_i} is {allele_j} : Removing");
                            remove.Add((index_i, allele_i));
                            colIndex[allele_i] = allele_j;
                        }
                        else
                        {
                            if (verbose)
                                Console.Error.WriteLine($"\t\t {allele_j} is {allele_i} : Removing");
                            remove.Add((index_j, allele_j));
                            colIndex[allele_j] = allele_i;
                        }
                        break;
                    }

                    if (seq_i_strip.Length < seq_j_strip.Length && seq_j_strip.Contains(seq_i_strip))
                    {
                        if (verbose)
                            Console.Error.WriteLine($"\t\t Collapsing {allele_i} into {allele_j}");
                        remove.Add((index_i, allele_i));
                        colIndex[allele_i] = allele_j;
                        break;
                    }
                }
            }

            // Remove from end to preserve indices
            remove.Sort((a, b) => b.Item1.CompareTo(a.Item1));

            foreach (var (idx, name) in remove)
            {
                seqs.RemoveAt(idx);
                index.Remove(name);
            }

            // Update indices
            var updatedIndex = new Dictionary<string, int>();
            foreach (var kvp in index)
            {
                int ind = kvp.Value;
                foreach (var (removedIdx, _) in remove)
                {
                    if (ind > removedIdx) ind--;
                }
                updatedIndex[kvp.Key] = ind;
            }

            return (updatedIndex, seqs, colIndex);
        }

        // Database operations
        public static void DownloadGenomeAndIndex(string destination)
        {
            if (!Directory.Exists(destination))
                throw new DirectoryNotFoundException($"Destination does not exist: {destination}");

            string[] requiredFiles = {
                $"{destination}/grch38",
                $"{destination}/genome.fa",
                $"{destination}/genome.fa.fai"
            };

            if (requiredFiles.All(File.Exists))
                return;

            Console.Error.WriteLine("Downloading GRCh38 genome and HISAT2 indexes...");

            var commands = new[]
            {
                "wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grch38.tar.gz",
                "tar xvzf grch38.tar.gz",
                "rm grch38.tar.gz",
                "hisat2-inspect grch38/genome > genome.fa",
                "samtools faidx genome.fa",
                $"mv grch38 {destination}",
                $"mv genome.fa genome.fa.fai {destination}"
            };

            foreach (var cmd in commands)
            {
                Process.Start("bash", $"-c \"{cmd}\"")?.WaitForExit();
            }
        }

        public static void CloneHisatGenotypeDatabase(string destination)
        {
            string dbPath = Path.Combine(destination, "hisatgenotype_db");
            if (Directory.Exists(dbPath))
                return;

            Console.Error.WriteLine("Cloning hisatgenotype_db...");
            Process.Start("git", "clone https://github.com/DaehwanKimLab/hisatgenotype_db.git")?.WaitForExit();

            if (Directory.Exists("hisatgenotype_db"))
            {
                Directory.Move("hisatgenotype_db", dbPath);
            }
        }

        // Read genome sequence from FASTA file
        public static (Dictionary<string, string> chrDic, List<string> chrNames, List<string> chrFullNames)
            ReadGenome(string genomeFile)
        {
            var chrDic = new Dictionary<string, string>();
            var chrNames = new List<string>();
            var chrFullNames = new List<string>();

            string content = File.ReadAllText(genomeFile);
            var seqs = content.Trim('\n').Split('>').Skip(1).ToArray();

            foreach (var seqBlock in seqs)
            {
                int ix = seqBlock.IndexOf('\n');
                string chrFullName = seqBlock.Substring(0, ix);
                string sequence = seqBlock.Substring(ix).Replace("\n", "");
                string chrName = chrFullName.Split()[0];

                chrDic[chrName] = sequence;
                chrNames.Add(chrName);
                chrFullNames.Add(chrFullName);
            }

            return (chrDic, chrNames, chrFullNames);
        }

        // Write sequences to FASTA file
        public static void WriteFasta(string filename, Dictionary<string, string> sequences, bool addLen = true)
        {
            string filePath = Path.GetDirectoryName(filename);
            if (!string.IsNullOrEmpty(filePath) && !Directory.Exists(filePath))
            {
                Directory.CreateDirectory(filePath);
            }

            using (var writer = new StreamWriter(filename, append: true))
            {
                foreach (var kvp in sequences)
                {
                    if (addLen)
                        writer.WriteLine($">{kvp.Key} {kvp.Value.Length}bp");
                    else
                        writer.WriteLine($">{kvp.Key}");

                    for (int i = 0; i < kvp.Value.Length; i += 60)
                    {
                        int len = Math.Min(60, kvp.Value.Length - i);
                        writer.WriteLine(kvp.Value.Substring(i, len));
                    }
                }
            }
        }

        // Read variants from .snp file
        public static (Dictionary<string, Dictionary<string, (string type, int pos, string data)>> varData,
                      Dictionary<string, List<(int pos, string varId)>> varList)
            ReadVariants(string fname, bool genes = false)
        {
            var varData = new Dictionary<string, Dictionary<string, (string, int, string)>>();
            var varList = new Dictionary<string, List<(int, string)>>();

            foreach (var line in File.ReadLines(fname))
            {
                var parts = line.Trim().Split('\t');
                string varId = parts[0];
                string varType = parts[1];
                string name = parts[2];
                int pos = int.Parse(parts[3]);
                string varStr = parts[4];

                string gene = genes ? name.Split('*')[0] : name;

                if (!varData.ContainsKey(gene))
                {
                    varData[gene] = new Dictionary<string, (string, int, string)>();
                    varList[gene] = new List<(int, string)>();
                }

                if (genes)
                {
                    varData[gene][varId] = (varType, pos, varStr);
                    varList[gene].Add((pos, varId));
                }
                else
                {
                    // For non-gene mode, we still only track position and varId in the list
                    // The full data is in varData
                    if (!varData[gene].ContainsKey(varId))
                    {
                        varData[gene][varId] = (varType, pos, varStr);
                    }
                    varList[gene].Add((pos, varId));
                }
            }

            foreach (var gene in varList.Keys.ToList())
            {
                varList[gene] = varList[gene].OrderBy(x => x.Item1).ToList();
            }

            return (varData, varList);
        }

        // Read haplotypes from .haplotype file
        public static Dictionary<string, List<(int left, int right, List<string> vars)>>
            ReadHaplotypes(string fname)
        {
            var alleleHaplotypes = new Dictionary<string, List<(int, int, List<string>)>>();

            foreach (var line in File.ReadLines(fname))
            {
                var parts = line.Trim().Split();
                string haplotypeId = parts[0];
                string alleleName = parts[1];
                int left = int.Parse(parts[2]);
                int right = int.Parse(parts[3]);
                var vars = parts[4].Split(',').ToList();

                if (!alleleHaplotypes.ContainsKey(alleleName))
                {
                    alleleHaplotypes[alleleName] = new List<(int, int, List<string>)>();
                }

                alleleHaplotypes[alleleName].Add((left, right, vars));
            }

            return alleleHaplotypes;
        }

        // Read links from .link file
        public static object ReadLinks(string fname, bool asList = false)
        {
            if (asList)
            {
                var linksList = new List<(string varId, List<string> alleleNames)>();
                foreach (var line in File.ReadLines(fname))
                {
                    var parts = line.Trim().Replace(" ", "\t").Split('\t');
                    string varId = parts[0];
                    var alleleNames = parts.Skip(1).ToList();
                    linksList.Add((varId, alleleNames));
                }
                return linksList;
            }
            else
            {
                var linksDict = new Dictionary<string, List<string>>();
                foreach (var line in File.ReadLines(fname))
                {
                    var parts = line.Trim().Replace(" ", "\t").Split('\t');
                    string varId = parts[0];
                    var alleleNames = parts.Skip(1).ToList();
                    linksDict[varId] = alleleNames;
                }
                return linksDict;
            }
        }

        // Binary search the variant list to find the lower position
        public static int LowerBound(List<(int pos, string varId)> varList, int pos)
        {
            int low = 0, high = varList.Count;

            while (low < high)
            {
                int m = (low + high) / 2;
                int mPos = varList[m].pos;

                if (mPos < pos)
                {
                    low = m + 1;
                }
                else if (mPos > pos)
                {
                    high = m;
                }
                else
                {
                    // Found exact match, find leftmost
                    while (m > 0 && varList[m - 1].pos == pos)
                    {
                        m--;
                    }
                    return m;
                }
            }

            return low;
        }

        // Read .locus file for gene coordinates
        public static (Dictionary<string, string> refGenes,
                      Dictionary<string, (string geneName, string chrom, int left, int right,
                                         List<(int, int)> exons, List<(int, int)> primaryExons)> refGeneLoci)
            ReadLocus(string fname, bool isGenome, string target)
        {
            var refGenes = new Dictionary<string, string>();
            var refGeneLoci = new Dictionary<string, (string, string, int, int, List<(int, int)>, List<(int, int)>)>();

            foreach (var line in File.ReadLines(fname))
            {
                var locus = line.Trim().Split();
                string gene, geneName, chrom, exonStr, strand;
                int left, right;

                if (isGenome)
                {
                    gene = locus[0];
                    geneName = locus[1];
                    chrom = locus[2];
                    left = int.Parse(locus[3]);
                    right = int.Parse(locus[4]);
                    exonStr = locus[5];
                    strand = locus[6];

                    if (gene.ToLower() != target)
                        continue;
                }
                else
                {
                    geneName = locus[0];
                    chrom = locus[1];
                    left = int.Parse(locus[2]);
                    right = int.Parse(locus[3]);
                    exonStr = locus[5];
                    strand = locus[6];
                }

                string geneGene = geneName.Split('*')[0];
                refGenes[geneGene] = geneName;

                var exons = new List<(int, int)>();
                var primaryExons = new List<(int, int)>();

                foreach (var exon in exonStr.Split(','))
                {
                    bool primary = exon.EndsWith('p');
                    string exonClean = primary ? exon.Substring(0, exon.Length - 1) : exon;

                    var exonParts = exonClean.Split('-');
                    int exonLeft = int.Parse(exonParts[0]);
                    int exonRight = int.Parse(exonParts[1]);

                    exons.Add((exonLeft, exonRight));
                    if (primary)
                        primaryExons.Add((exonLeft, exonRight));
                }

                refGeneLoci[geneGene] = (geneName, chrom, left, right, exons, primaryExons);
            }

            return (refGenes, refGeneLoci);
        }

        // Read allele sequences from FASTA file
        public static Dictionary<string, string> ReadAlleleSeq(string fname, bool genes = false)
        {
            var dic = new Dictionary<string, string>();

            string content = File.ReadAllText(fname);
            var seqs = content.Trim('\n').Split('>').Skip(1).ToArray();

            foreach (var seqBlock in seqs)
            {
                int ix = seqBlock.IndexOf('\n');
                string seqName = seqBlock.Substring(0, ix);
                string sequence = seqBlock.Substring(ix).Replace("\n", "");

                if (genes)
                {
                    string gene = seqName.Split('*')[0];
                    if (!dic.ContainsKey(gene))
                    {
                        dic[gene] = sequence;
                    }
                }
                else
                {
                    dic[seqName] = sequence;
                }
            }

            return dic;
        }
    }
}
