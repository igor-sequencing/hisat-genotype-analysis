// ---------------------------------------------------------------------------
// Copyright 2015, Daehwan Kim <infphilo@gmail.com>
// Ported to C# 2025
//
// This file is part of HISAT 2. This script is a wrapper for the functions
// needed to build HISATgenotype indices
//
// HISAT 2 is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// HISAT 2 is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with HISAT 2.  If not, see <http://www.gnu.org/licenses/>.
// ---------------------------------------------------------------------------

using System;
using System.Collections.Generic;
using System.CommandLine;
using System.CommandLine.Invocation;
using System.IO;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;

namespace HisatGenotype.Tools
{
    /// <summary>
    /// Extract variants from multiple sequence alignments
    /// </summary>
    public class HisatGenotypeExtractVars
    {
        private static readonly SemaphoreSlim _lock = new SemaphoreSlim(1, 1);

        public static async Task<int> Main(string[] args)
        {
            var rootCommand = new RootCommand("Extract variants from multiple sequence alignments");

            // Add Arguments - Database options
            var baseNameOption = new Option<string>(
                aliases: new[] { "--base", "--base-fname" },
                description: "Base file name for index, variants, haplotypes, etc. (e.g. hla, rbg, codis)",
                getDefaultValue: () => "");

            var locusListOption = new Option<string>(
                aliases: new[] { "--locus-list" },
                description: "A comma-separated list of gene names (default: empty, all genes)",
                getDefaultValue: () => "");

            var indexDirOption = new Option<string>(
                aliases: new[] { "-z", "--index-dir" },
                description: "Set location to use for indices",
                getDefaultValue: () => GetDefaultIndexDirectory());

            // Variable gap options
            var interGapOption = new Option<int>(
                aliases: new[] { "--inter-gap" },
                description: "Maximum distance for variants to be in the same haplotype",
                getDefaultValue: () => 30);

            var intraGapOption = new Option<int>(
                aliases: new[] { "--intra-gap" },
                description: "Break a haplotype into several haplotypes",
                getDefaultValue: () => 50);

            // Extract vars options
            var wholeHaplotypeOption = new Option<bool>(
                aliases: new[] { "--whole-haplotype" },
                description: "Include partial alleles (e.g. A_nuc.fasta)",
                getDefaultValue: () => false);

            var minVarFreqOption = new Option<double>(
                aliases: new[] { "--min-var-freq" },
                description: "Exclude variants whose freq is below than this value in percentage",
                getDefaultValue: () => 0.0);

            var extSeqOption = new Option<int>(
                aliases: new[] { "--ext-seq" },
                description: "Length of extra sequences flanking backbone sequences",
                getDefaultValue: () => 0);

            var leftshiftOption = new Option<bool>(
                aliases: new[] { "--leftshift" },
                description: "Shift deletions to the leftmost",
                getDefaultValue: () => false);

            var partialOption = new Option<bool>(
                aliases: new[] { "--no-partial" },
                description: "Include partial alleles (e.g. A_nuc.fasta)",
                getDefaultValue: () => true);

            // Common options
            var threadsOption = new Option<int>(
                aliases: new[] { "-p", "--threads" },
                description: "Number of threads",
                getDefaultValue: () => 1);

            var verboseOption = new Option<bool>(
                aliases: new[] { "-v", "--verbose" },
                description: "Print statistics to stderr",
                getDefaultValue: () => false);

            // Add all options to command
            rootCommand.AddOption(baseNameOption);
            rootCommand.AddOption(locusListOption);
            rootCommand.AddOption(indexDirOption);
            rootCommand.AddOption(interGapOption);
            rootCommand.AddOption(intraGapOption);
            rootCommand.AddOption(wholeHaplotypeOption);
            rootCommand.AddOption(minVarFreqOption);
            rootCommand.AddOption(extSeqOption);
            rootCommand.AddOption(leftshiftOption);
            rootCommand.AddOption(partialOption);
            rootCommand.AddOption(threadsOption);
            rootCommand.AddOption(verboseOption);

            rootCommand.SetHandler(async (context) =>
            {
                var baseName = context.ParseResult.GetValueForOption(baseNameOption);
                var locusList = context.ParseResult.GetValueForOption(locusListOption);
                var indexDir = context.ParseResult.GetValueForOption(indexDirOption);
                var interGap = context.ParseResult.GetValueForOption(interGapOption);
                var intraGap = context.ParseResult.GetValueForOption(intraGapOption);
                var wholeHaplotype = context.ParseResult.GetValueForOption(wholeHaplotypeOption);
                var minVarFreq = context.ParseResult.GetValueForOption(minVarFreqOption);
                var extSeq = context.ParseResult.GetValueForOption(extSeqOption);
                var leftshift = context.ParseResult.GetValueForOption(leftshiftOption);
                var partial = context.ParseResult.GetValueForOption(partialOption);
                var threads = context.ParseResult.GetValueForOption(threadsOption);
                var verbose = context.ParseResult.GetValueForOption(verboseOption);

                try
                {
                    await ExecuteExtractVars(
                        baseName,
                        locusList,
                        indexDir,
                        interGap,
                        intraGap,
                        wholeHaplotype,
                        minVarFreq,
                        extSeq,
                        leftshift,
                        partial,
                        threads,
                        verbose);

                    context.ExitCode = 0;
                }
                catch (Exception ex)
                {
                    Console.Error.WriteLine($"Error: {ex.Message}");
                    Console.Error.WriteLine(ex.StackTrace);
                    context.ExitCode = 1;
                }
            });

            return await rootCommand.InvokeAsync(args);
        }

        private static string GetDefaultIndexDirectory()
        {
            var currentDir = Directory.GetCurrentDirectory();
            var linkFile = Path.Combine(currentDir, "hg_ix.link");

            if (File.Exists(linkFile))
            {
                return File.ReadAllText(linkFile).Trim();
            }

            return Path.Combine(currentDir, "indicies");
        }

        private static async Task ExecuteExtractVars(
            string baseName,
            string locusList,
            string indexDir,
            int interGap,
            int intraGap,
            bool wholeHaplotype,
            double minVarFreq,
            int extSeq,
            bool leftshift,
            bool partial,
            int threads,
            bool verbose)
        {
            // Validate parameters
            if (interGap > intraGap)
            {
                throw new ArgumentException(
                    $"--inter-gap ({interGap}) must be smaller than --intra-gap ({intraGap})");
            }

            // Parse locus list
            var loci = string.IsNullOrEmpty(locusList)
                ? new List<string>()
                : locusList.Split(',').Select(s => s.Trim()).ToList();

            // Clone hisatgenotype database from git if needed
            var dbPath = Path.Combine(indexDir, "hisatgenotype_db");
            if (!Directory.Exists(dbPath))
            {
                CloneHisatGenotypeDatabase(indexDir);
            }

            // Get base file names
            List<string> baseNames;
            if (string.IsNullOrEmpty(baseName))
            {
                baseNames = new List<string>();
                foreach (var dir in Directory.GetDirectories(dbPath))
                {
                    var dirName = Path.GetFileName(dir);
                    if (dirName != ".git" && dirName != "README.md")
                    {
                        baseNames.Add(dirName.ToLower());
                    }
                }
            }
            else
            {
                baseNames = baseName.ToLower().Split(',').Select(s => s.Trim()).ToList();
            }

            // Process each database in parallel
            var options = new ParallelOptions
            {
                MaxDegreeOfParallelism = threads
            };

            var exceptions = new List<Exception>();

            await Task.Run(() =>
            {
                Parallel.ForEach(baseNames, options, (baseFile) =>
                {
                    try
                    {
                        string baseDirName = "";
                        string baseFileName = baseFile;

                        if (baseFile.Contains('/'))
                        {
                            var parts = baseFile.Split('/');
                            baseFileName = parts[^1];
                            baseDirName = string.Join('/', parts[..^1]);
                        }

                        var targetDir = string.IsNullOrEmpty(baseDirName) ? indexDir : baseDirName;

                        ExtractVarsForBase(
                            baseFileName,
                            targetDir,
                            loci,
                            interGap,
                            intraGap,
                            wholeHaplotype,
                            minVarFreq,
                            extSeq,
                            leftshift,
                            partial,
                            verbose);
                    }
                    catch (Exception ex)
                    {
                        lock (exceptions)
                        {
                            exceptions.Add(ex);
                        }
                    }
                });
            });

            // Check for errors
            if (exceptions.Any())
            {
                throw new AggregateException("Errors occurred during variant extraction:", exceptions);
            }
        }

        private static void ExtractVarsForBase(
            string baseName,
            string indexDir,
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
            // NOTE: This method would call into the actual variant extraction logic
            // from hisatgenotype_typing_process module. In a full C# port, that module
            // would also need to be translated to C#.

            // For now, this is a wrapper that would call the Python module via
            // Process, or ideally call a native C# implementation.

            if (verbose)
            {
                Console.Error.WriteLine($"Processing base: {baseName}");
                Console.Error.WriteLine($"  Index directory: {indexDir}");
                Console.Error.WriteLine($"  Locus list: {string.Join(",", locusList)}");
                Console.Error.WriteLine($"  Inter-gap: {interGap}");
                Console.Error.WriteLine($"  Intra-gap: {intraGap}");
                Console.Error.WriteLine($"  Whole haplotype: {wholeHaplotype}");
                Console.Error.WriteLine($"  Min variant frequency: {minVarFreq}");
                Console.Error.WriteLine($"  Extended sequence length: {extSeqLen}");
                Console.Error.WriteLine($"  Leftshift: {leftshift}");
                Console.Error.WriteLine($"  Include partial: {partial}");
            }

            // TODO: Call the actual extract_vars implementation
            // This requires porting hisatgenotype_typing_process.extract_vars to C#

            // Example placeholder call:
            HisatGenotypeTypingProcess.ExtractVars(
                baseName,
                indexDir,
                locusList,
                interGap,
                intraGap,
                wholeHaplotype,
                minVarFreq,
                extSeqLen,
                leftshift,
                partial,
                verbose);
        }

        private static void CloneHisatGenotypeDatabase(string indexDir)
        {
            // TODO: Implement git clone functionality
            // This would use LibGit2Sharp or System.Diagnostics.Process to run git

            Console.Error.WriteLine($"Cloning hisatgenotype_db to {indexDir}...");

            var startInfo = new System.Diagnostics.ProcessStartInfo
            {
                FileName = "git",
                Arguments = "clone https://github.com/DaehwanKimLab/hisatgenotype_db.git",
                WorkingDirectory = indexDir,
                UseShellExecute = false,
                RedirectStandardOutput = true,
                RedirectStandardError = true
            };

            using var process = System.Diagnostics.Process.Start(startInfo);
            process?.WaitForExit();

            if (process?.ExitCode != 0)
            {
                throw new Exception($"Failed to clone hisatgenotype_db. Exit code: {process?.ExitCode}");
            }
        }
    }

    /// <summary>
    /// Placeholder class for the typing process module
    /// This would need to be fully ported from the Python module
    /// </summary>
    public static class HisatGenotypeTypingProcess
    {
        public static void ExtractVars(
            string baseName,
            string indexDir,
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
            // NOTE: This is a placeholder for the actual extract_vars implementation
            // The full Python module hisatgenotype_typing_process.py would need to be
            // ported to C# for complete functionality.

            throw new NotImplementedException(
                "The extract_vars function from hisatgenotype_typing_process module " +
                "needs to be ported to C#. This requires translating the Python module " +
                "hisatgenotype_modules/hisatgenotype_typing_process.py to C#.");
        }
    }
}
