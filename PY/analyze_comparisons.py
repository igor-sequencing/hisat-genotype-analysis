#!/usr/bin/env python3
"""
Analyze comparison results across all samples and generate summary statistics.
"""

import re
from pathlib import Path
from collections import defaultdict
from bs4 import BeautifulSoup

def parse_html_comparison(html_file):
    """Parse a comparison HTML file and extract allele calls."""
    with open(html_file, 'r') as f:
        soup = BeautifulSoup(f.read(), 'html.parser')

    sample_name = soup.find('strong', string=lambda x: x and 'strong' in str(soup.find('strong')))
    if not sample_name:
        # Try alternative method
        sample_info = soup.find('div', class_='sample-info')
        if sample_info:
            sample_name = sample_info.find('strong').text

    data = {
        'sample': Path(html_file).stem.replace('_comparison', ''),
        'genes': defaultdict(dict)
    }

    # Parse table rows
    rows = soup.find_all('tr')
    for row in rows[1:]:  # Skip header
        cells = row.find_all('td')
        if len(cells) < 6:
            continue

        gene = cells[0].text.strip()

        # OptiType
        optitype_alleles = [span.text.strip() for span in cells[1].find_all('span', class_='allele')]
        data['genes'][gene]['optitype'] = optitype_alleles if optitype_alleles else []

        # HLA-LA
        hlala_alleles = [span.text.strip() for span in cells[2].find_all('span', class_='allele')]
        data['genes'][gene]['hlala'] = hlala_alleles if hlala_alleles else []

        # HLA-HD
        hlahd_alleles = [span.text.strip() for span in cells[3].find_all('span', class_='allele')]
        data['genes'][gene]['hlahd'] = hlahd_alleles if hlahd_alleles else []

        # HISAT-genotype
        hisat_alleles = []
        for span in cells[4].find_all('span', class_='allele'):
            allele_text = span.text.strip()
            # Extract just the allele name before (rank:
            allele = allele_text.split(' (rank:')[0].strip() if '(rank:' in allele_text else allele_text
            hisat_alleles.append(allele)
        data['genes'][gene]['hisat'] = hisat_alleles if hisat_alleles else []

        # Dragen
        dragen_alleles = [span.text.strip() for span in cells[5].find_all('span', class_='allele')]
        data['genes'][gene]['dragen'] = dragen_alleles if dragen_alleles else []

    return data

def calculate_statistics(all_samples):
    """Calculate statistics across all samples."""
    stats = {
        'total_samples': len(all_samples),
        'genes_analyzed': set(),
        'method_coverage': defaultdict(lambda: defaultdict(int)),  # method -> gene -> count
        'method_allele_counts': defaultdict(int),  # method -> total alleles called
        'concordance': defaultdict(lambda: defaultdict(int)),  # gene -> method pairs -> matching count
        'avg_alleles_per_gene': defaultdict(lambda: defaultdict(list)),  # method -> gene -> [counts]
    }

    # Collect all genes
    for sample in all_samples:
        stats['genes_analyzed'].update(sample['genes'].keys())

    stats['genes_analyzed'] = sorted(stats['genes_analyzed'])

    # Calculate coverage and allele counts
    for sample in all_samples:
        for gene in stats['genes_analyzed']:
            if gene in sample['genes']:
                for method in ['optitype', 'hlala', 'hlahd', 'hisat', 'dragen']:
                    alleles = sample['genes'][gene].get(method, [])
                    if alleles:
                        stats['method_coverage'][method][gene] += 1
                        stats['method_allele_counts'][method] += len(alleles)
                        stats['avg_alleles_per_gene'][method][gene].append(len(alleles))

    # Calculate concordance between methods (exact allele matching)
    for sample in all_samples:
        for gene in sample['genes']:
            methods = sample['genes'][gene]

            # Compare each pair of methods
            method_pairs = [
                ('optitype', 'hlala'),
                ('optitype', 'hlahd'),
                ('optitype', 'hisat'),
                ('hlala', 'hlahd'),
                ('hlala', 'hisat'),
                ('hlahd', 'hisat')
            ]

            for m1, m2 in method_pairs:
                alleles1 = set(methods.get(m1, []))
                alleles2 = set(methods.get(m2, []))

                if alleles1 and alleles2:
                    # Check if there's any overlap
                    if alleles1 & alleles2:
                        stats['concordance'][gene][f'{m1}_vs_{m2}'] += 1

    return stats

def generate_summary_html(stats, all_samples, output_file):
    """Generate HTML summary report."""

    html = """<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>HLA Typing Methods Comparison Summary</title>
    <style>
        body {
            font-family: Arial, sans-serif;
            margin: 20px;
            background-color: #f5f5f5;
        }
        h1, h2 {
            color: #333;
        }
        .container {
            background-color: white;
            padding: 20px;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            margin-bottom: 20px;
        }
        table {
            border-collapse: collapse;
            width: 100%;
            margin: 20px 0;
        }
        th, td {
            border: 1px solid #ddd;
            padding: 12px;
            text-align: left;
        }
        th {
            background-color: #4CAF50;
            color: white;
            font-weight: bold;
        }
        tr:nth-child(even) {
            background-color: #f9f9f9;
        }
        tr:hover {
            background-color: #f1f1f1;
        }
        .stat-box {
            display: inline-block;
            margin: 10px;
            padding: 15px;
            background-color: #e3f2fd;
            border-radius: 5px;
            min-width: 200px;
        }
        .stat-label {
            font-size: 14px;
            color: #666;
        }
        .stat-value {
            font-size: 24px;
            font-weight: bold;
            color: #1976d2;
        }
        .method-optitype { background-color: #f3e5f5; }
        .method-hlala { background-color: #e8f5e9; }
        .method-hlahd { background-color: #fff3e0; }
        .method-hisat { background-color: #e3f2fd; }
        .method-dragen { background-color: #fce4ec; }
    </style>
</head>
<body>
    <h1>HLA Typing Methods Comparison Summary</h1>
"""

    # Overall statistics
    html += f"""
    <div class="container">
        <h2>Overall Statistics</h2>
        <div class="stat-box">
            <div class="stat-label">Total Samples</div>
            <div class="stat-value">{stats['total_samples']}</div>
        </div>
        <div class="stat-box">
            <div class="stat-label">Genes Analyzed</div>
            <div class="stat-value">{len(stats['genes_analyzed'])}</div>
        </div>
        <div class="stat-box">
            <div class="stat-label">Total Comparisons</div>
            <div class="stat-value">{stats['total_samples'] * len(stats['genes_analyzed'])}</div>
        </div>
    </div>
"""

    # Method coverage summary
    html += """
    <div class="container">
        <h2>Method Coverage by Gene</h2>
        <p>Number of samples where each method successfully called alleles for each gene</p>
        <table>
            <thead>
                <tr>
                    <th>Gene</th>
                    <th class="method-optitype">OptiType</th>
                    <th class="method-hlala">HLA-LA</th>
                    <th class="method-hlahd">HLA-HD</th>
                    <th class="method-hisat">HISAT-genotype</th>
                    <th class="method-dragen">Dragen 4.4</th>
                </tr>
            </thead>
            <tbody>
"""

    for gene in stats['genes_analyzed']:
        html += f"                <tr>\n"
        html += f"                    <td><strong>{gene}</strong></td>\n"
        html += f"                    <td>{stats['method_coverage']['optitype'][gene]} / {stats['total_samples']}</td>\n"
        html += f"                    <td>{stats['method_coverage']['hlala'][gene]} / {stats['total_samples']}</td>\n"
        html += f"                    <td>{stats['method_coverage']['hlahd'][gene]} / {stats['total_samples']}</td>\n"
        html += f"                    <td>{stats['method_coverage']['hisat'][gene]} / {stats['total_samples']}</td>\n"
        html += f"                    <td>{stats['method_coverage']['dragen'][gene]} / {stats['total_samples']}</td>\n"
        html += f"                </tr>\n"

    html += """            </tbody>
        </table>
    </div>
"""

    # Total alleles called
    html += """
    <div class="container">
        <h2>Total Alleles Called</h2>
        <table>
            <thead>
                <tr>
                    <th>Method</th>
                    <th>Total Alleles</th>
                    <th>Average per Sample</th>
                </tr>
            </thead>
            <tbody>
"""

    for method in ['optitype', 'hlala', 'hlahd', 'hisat', 'dragen']:
        total = stats['method_allele_counts'][method]
        avg = total / stats['total_samples'] if stats['total_samples'] > 0 else 0
        method_name = {
            'optitype': 'OptiType',
            'hlala': 'HLA-LA',
            'hlahd': 'HLA-HD',
            'hisat': 'HISAT-genotype',
            'dragen': 'Dragen 4.4'
        }[method]
        html += f"                <tr>\n"
        html += f"                    <td><strong>{method_name}</strong></td>\n"
        html += f"                    <td>{total}</td>\n"
        html += f"                    <td>{avg:.1f}</td>\n"
        html += f"                </tr>\n"

    html += """            </tbody>
        </table>
    </div>
"""

    # Average alleles per gene
    html += """
    <div class="container">
        <h2>Average Alleles per Gene (when called)</h2>
        <table>
            <thead>
                <tr>
                    <th>Gene</th>
                    <th class="method-optitype">OptiType</th>
                    <th class="method-hlala">HLA-LA</th>
                    <th class="method-hlahd">HLA-HD</th>
                    <th class="method-hisat">HISAT-genotype</th>
                    <th class="method-dragen">Dragen 4.4</th>
                </tr>
            </thead>
            <tbody>
"""

    for gene in stats['genes_analyzed']:
        html += f"                <tr>\n"
        html += f"                    <td><strong>{gene}</strong></td>\n"
        for method in ['optitype', 'hlala', 'hlahd', 'hisat', 'dragen']:
            counts = stats['avg_alleles_per_gene'][method][gene]
            avg = sum(counts) / len(counts) if counts else 0
            html += f"                    <td>{avg:.2f}</td>\n"
        html += f"                </tr>\n"

    html += """            </tbody>
        </table>
    </div>
"""

    html += """
</body>
</html>
"""

    with open(output_file, 'w') as f:
        f.write(html)

    print(f"Generated summary: {output_file}")

def main():
    comparison_dir = Path('comparison_results')
    output_file = 'hla_methods_summary.html'

    if not comparison_dir.exists():
        print(f"Error: {comparison_dir} not found")
        return

    # Parse all HTML files
    html_files = list(comparison_dir.glob('*_comparison.html'))

    if not html_files:
        print(f"No comparison HTML files found in {comparison_dir}")
        return

    print(f"Analyzing {len(html_files)} comparison files...")

    all_samples = []
    for html_file in sorted(html_files):
        print(f"  Parsing {html_file.name}...")
        data = parse_html_comparison(html_file)
        all_samples.append(data)

    # Calculate statistics
    print("\nCalculating statistics...")
    stats = calculate_statistics(all_samples)

    # Generate summary
    print("Generating summary report...")
    generate_summary_html(stats, all_samples, output_file)

    # Print text summary
    print("\n" + "="*60)
    print("SUMMARY STATISTICS")
    print("="*60)
    print(f"\nTotal samples analyzed: {stats['total_samples']}")
    print(f"Total genes analyzed: {len(stats['genes_analyzed'])}")
    print(f"\nGenes: {', '.join(stats['genes_analyzed'])}")

    print("\n" + "-"*60)
    print("Method Coverage Summary:")
    print("-"*60)
    for method in ['optitype', 'hlala', 'hlahd', 'hisat', 'dragen']:
        method_name = {
            'optitype': 'OptiType',
            'hlala': 'HLA-LA',
            'hlahd': 'HLA-HD',
            'hisat': 'HISAT-genotype',
            'dragen': 'Dragen 4.4'
        }[method]
        total = stats['method_allele_counts'][method]
        genes_covered = len([g for g in stats['genes_analyzed'] if stats['method_coverage'][method][g] > 0])
        print(f"{method_name:20s}: {total:4d} total alleles, {genes_covered:2d} genes covered")

    print(f"\n✓ Done! Open {output_file} to view the full report")

if __name__ == '__main__':
    main()
