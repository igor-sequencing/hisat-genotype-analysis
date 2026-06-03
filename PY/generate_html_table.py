#!/usr/bin/env python3
"""
Convert HISAT-genotype results to HTML table.
Rows: genes, Columns: sample directories, Values: ranked alleles with abundance
"""

import os
import re
from pathlib import Path
from collections import defaultdict

def parse_report_file(report_path):
    """Parse a single report file and extract ranked results by gene."""
    results = {}

    with open(report_path, 'r') as f:
        content = f.read()

    # Find all ranked sections using regex
    # Pattern: "1 ranked GENE*ALLELE (abundance: XX.XX%)"
    pattern = r'(\d+) ranked ([A-Z0-9]+)\*([^\s]+) \(abundance: ([\d.]+)%\)'
    matches = re.findall(pattern, content)

    current_gene = None
    for rank, gene, allele, abundance in matches:
        if gene not in results:
            results[gene] = []
        results[gene].append({
            'rank': int(rank),
            'allele': allele,
            'abundance': float(abundance)
        })

    return results

def collect_all_results(output_dir):
    """Collect results from all directories in the output folder."""
    all_data = {}
    samples = []

    # Get all subdirectories
    output_path = Path(output_dir)
    if not output_path.exists():
        print(f"Error: Directory {output_dir} does not exist")
        return None, None, None

    for sample_dir in sorted(output_path.iterdir()):
        if not sample_dir.is_dir():
            continue

        sample_name = sample_dir.name

        # Find report file in this directory
        report_files = list(sample_dir.glob('*.report'))
        if not report_files:
            print(f"Warning: No report file found in {sample_name}")
            continue

        report_file = report_files[0]
        results = parse_report_file(report_file)

        if results:
            samples.append(sample_name)
            all_data[sample_name] = results

    # Get all unique genes across all samples
    all_genes = set()
    for sample_results in all_data.values():
        all_genes.update(sample_results.keys())

    return all_data, sorted(samples), sorted(all_genes)

def generate_html_table(all_data, samples, genes, output_file):
    """Generate HTML table with genes as rows and samples as columns."""

    html = """<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>HISAT-genotype Results</title>
    <style>
        body {
            font-family: Arial, sans-serif;
            margin: 20px;
            background-color: #f5f5f5;
        }
        h1 {
            color: #333;
            text-align: center;
        }
        .table-container {
            overflow-x: auto;
            background-color: white;
            padding: 20px;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }
        table {
            border-collapse: collapse;
            width: 100%;
            min-width: 800px;
        }
        th, td {
            border: 1px solid #ddd;
            padding: 8px;
            text-align: left;
            font-size: 12px;
        }
        th {
            background-color: #4CAF50;
            color: white;
            position: sticky;
            top: 0;
            z-index: 10;
        }
        th.gene-col {
            background-color: #2196F3;
            min-width: 100px;
        }
        tr:nth-child(even) {
            background-color: #f2f2f2;
        }
        tr:hover {
            background-color: #ddd;
        }
        .allele {
            display: block;
            margin: 2px 0;
        }
        .rank-1 { color: #d32f2f; font-weight: bold; }
        .rank-2 { color: #f57c00; font-weight: bold; }
        .rank-3 { color: #7b1fa2; }
        .rank-4 { color: #1976d2; }
        .rank-5 { color: #388e3c; }
    </style>
</head>
<body>
    <h1>HISAT-genotype Results Summary</h1>
    <div class="table-container">
        <table>
            <thead>
                <tr>
                    <th class="gene-col">Gene</th>
"""

    # Add sample column headers
    for sample in samples:
        html += f"                    <th>{sample}</th>\n"

    html += """                </tr>
            </thead>
            <tbody>
"""

    # Add rows for each gene
    for gene in genes:
        html += f"                <tr>\n                    <td><strong>{gene}</strong></td>\n"

        for sample in samples:
            html += "                    <td>"

            if sample in all_data and gene in all_data[sample]:
                ranked_alleles = all_data[sample][gene]
                # Sort by rank
                ranked_alleles.sort(key=lambda x: x['rank'])

                # Display each ranked allele
                for allele_info in ranked_alleles:
                    rank = allele_info['rank']
                    allele = allele_info['allele']
                    abundance = allele_info['abundance']
                    html += f'<span class="allele rank-{rank}">{rank}. {allele} ({abundance}%)</span>'
            else:
                html += "-"

            html += "</td>\n"

        html += "                </tr>\n"

    html += """            </tbody>
        </table>
    </div>
</body>
</html>
"""

    with open(output_file, 'w') as f:
        f.write(html)

    print(f"HTML table generated: {output_file}")

def main():
    output_dir = 'hisat.output'
    output_file = 'hisat_results_table.html'

    print(f"Parsing results from {output_dir}...")
    all_data, samples, genes = collect_all_results(output_dir)

    if not all_data:
        print("No data found!")
        return

    print(f"Found {len(samples)} samples and {len(genes)} genes")
    print(f"Generating HTML table...")

    generate_html_table(all_data, samples, genes, output_file)
    print(f"\nDone! Open {output_file} in your browser to view the results.")

if __name__ == '__main__':
    main()
