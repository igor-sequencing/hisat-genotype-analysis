#!/usr/bin/env python3
"""
Analyze concordance between HLA typing methods to assess reliability.
"""

import re
from pathlib import Path
from collections import defaultdict
from bs4 import BeautifulSoup

def normalize_allele(allele):
    """Normalize allele names to different resolution levels for comparison."""
    if not allele or allele == '-':
        return None

    # Remove HLA- prefix if present
    allele = allele.replace('HLA-', '')

    # Split by * to get gene and allele
    parts = allele.split('*')
    if len(parts) < 2:
        return None

    gene = parts[0]
    allele_code = parts[1]

    # Split allele code by : to get different resolution levels
    fields = allele_code.split(':')

    return {
        'full': allele,
        'gene': gene,
        '2field': f"{gene}*{':'.join(fields[:2])}" if len(fields) >= 2 else allele,
        '1field': f"{gene}*{fields[0]}" if len(fields) >= 1 else allele,
    }

def parse_html_comparison(html_file):
    """Parse a comparison HTML file and extract allele calls."""
    with open(html_file, 'r') as f:
        soup = BeautifulSoup(f.read(), 'html.parser')

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

def calculate_concordance(all_samples):
    """Calculate concordance statistics between methods."""

    stats = {
        'by_gene': defaultdict(lambda: {
            'total_samples': 0,
            'all_4_methods': 0,
            'any_3_methods': 0,
            'any_2_methods': 0,
            'only_1_method': 0,
            'concordance_2field': defaultdict(int),  # How many methods agree at 2-field
            'concordance_1field': defaultdict(int),  # How many methods agree at 1-field
            'full_agreement': 0,
        }),
        'overall': {
            'total_comparisons': 0,
            'all_4_agree_2field': 0,
            'any_3_agree_2field': 0,
            'any_2_agree_2field': 0,
            'no_agreement': 0,
        },
        'method_pairs': defaultdict(lambda: {'compared': 0, 'agree_2field': 0, 'agree_1field': 0}),
    }

    methods = ['optitype', 'hlala', 'hlahd', 'hisat', 'dragen']
    method_pairs = [
        ('optitype', 'hlala'),
        ('optitype', 'hlahd'),
        ('optitype', 'hisat'),
        ('optitype', 'dragen'),
        ('hlala', 'hlahd'),
        ('hlala', 'hisat'),
        ('hlala', 'dragen'),
        ('hlahd', 'hisat'),
        ('hlahd', 'dragen'),
        ('hisat', 'dragen')
    ]

    for sample in all_samples:
        for gene, gene_data in sample['genes'].items():
            # Count how many methods have data
            methods_with_data = [m for m in methods if gene_data.get(m)]
            num_methods = len(methods_with_data)

            if num_methods == 0:
                continue

            stats['by_gene'][gene]['total_samples'] += 1

            if num_methods == 5:
                stats['by_gene'][gene]['all_4_methods'] += 1  # Keep var name for now
            elif num_methods == 4:
                stats['by_gene'][gene]['all_4_methods'] += 1
            elif num_methods == 3:
                stats['by_gene'][gene]['any_3_methods'] += 1
            elif num_methods == 2:
                stats['by_gene'][gene]['any_2_methods'] += 1
            elif num_methods == 1:
                stats['by_gene'][gene]['only_1_method'] += 1

            # Normalize all alleles for comparison
            normalized = {}
            for method in methods:
                alleles = gene_data.get(method, [])
                normalized[method] = []
                for allele in alleles:
                    norm = normalize_allele(allele)
                    if norm:
                        normalized[method].append(norm)

            # Check concordance at different resolutions
            if num_methods >= 2:
                stats['overall']['total_comparisons'] += 1

                # Get all 2-field alleles from all methods
                all_2field = set()
                for method in methods_with_data:
                    for norm in normalized[method]:
                        all_2field.add(norm['2field'])

                # Count how many methods have each allele
                agreement_count = 0
                for allele_2field in all_2field:
                    methods_with_allele = 0
                    for method in methods_with_data:
                        if any(norm['2field'] == allele_2field for norm in normalized[method]):
                            methods_with_allele += 1
                    agreement_count = max(agreement_count, methods_with_allele)
                    stats['by_gene'][gene]['concordance_2field'][methods_with_allele] += 1

                # Overall stats
                if agreement_count >= 4:
                    stats['overall']['all_4_agree_2field'] += 1
                    stats['by_gene'][gene]['full_agreement'] += 1
                elif agreement_count >= 3:
                    stats['overall']['any_3_agree_2field'] += 1
                elif agreement_count == 2:
                    stats['overall']['any_2_agree_2field'] += 1
                else:
                    stats['overall']['no_agreement'] += 1

            # Pairwise comparisons
            for m1, m2 in method_pairs:
                if gene_data.get(m1) and gene_data.get(m2):
                    stats['method_pairs'][f'{m1}_vs_{m2}']['compared'] += 1

                    # Check for any overlap at 2-field level
                    alleles1_2f = set(norm['2field'] for norm in normalized[m1])
                    alleles2_2f = set(norm['2field'] for norm in normalized[m2])

                    if alleles1_2f & alleles2_2f:
                        stats['method_pairs'][f'{m1}_vs_{m2}']['agree_2field'] += 1

                    # Check at 1-field level
                    alleles1_1f = set(norm['1field'] for norm in normalized[m1])
                    alleles2_1f = set(norm['1field'] for norm in normalized[m2])

                    if alleles1_1f & alleles2_1f:
                        stats['method_pairs'][f'{m1}_vs_{m2}']['agree_1field'] += 1

    return stats

def generate_concordance_html(stats, output_file):
    """Generate HTML report for concordance analysis."""

    html = """<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>HLA Typing Methods Concordance Analysis</title>
    <style>
        body {
            font-family: Arial, sans-serif;
            margin: 20px;
            background-color: #f5f5f5;
        }
        h1, h2, h3 {
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
        .high-concordance { background-color: #c8e6c9; }
        .medium-concordance { background-color: #fff9c4; }
        .low-concordance { background-color: #ffccbc; }
        .highlight {
            background-color: #fff3e0;
            padding: 15px;
            border-left: 4px solid #ff9800;
            margin: 15px 0;
        }
    </style>
</head>
<body>
    <h1>HLA Typing Methods Concordance Analysis</h1>

    <div class="highlight">
        <strong>Key Question:</strong> How reliable is allele identification when using multiple methods together?<br>
        <strong>Analysis:</strong> This report shows agreement between OptiType, HLA-LA, HLA-HD, and HISAT-genotype at 2-field resolution (e.g., A*02:01)
    </div>
"""

    # Overall concordance summary
    total = stats['overall']['total_comparisons']
    if total > 0:
        all4_pct = (stats['overall']['all_4_agree_2field'] / total) * 100
        any3_pct = (stats['overall']['any_3_agree_2field'] / total) * 100
        any2_pct = (stats['overall']['any_2_agree_2field'] / total) * 100
        no_pct = (stats['overall']['no_agreement'] / total) * 100
    else:
        all4_pct = any3_pct = any2_pct = no_pct = 0

    html += f"""
    <div class="container">
        <h2>Overall Concordance Summary</h2>
        <p>Total gene-sample combinations analyzed: {total}</p>
        <div class="stat-box high-concordance">
            <div class="stat-label">All 4 Methods Agree</div>
            <div class="stat-value">{stats['overall']['all_4_agree_2field']} ({all4_pct:.1f}%)</div>
        </div>
        <div class="stat-box high-concordance">
            <div class="stat-label">3+ Methods Agree</div>
            <div class="stat-value">{stats['overall']['any_3_agree_2field']} ({any3_pct:.1f}%)</div>
        </div>
        <div class="stat-box medium-concordance">
            <div class="stat-label">2 Methods Agree</div>
            <div class="stat-value">{stats['overall']['any_2_agree_2field']} ({any2_pct:.1f}%)</div>
        </div>
        <div class="stat-box low-concordance">
            <div class="stat-label">No Agreement</div>
            <div class="stat-value">{stats['overall']['no_agreement']} ({no_pct:.1f}%)</div>
        </div>
    </div>
"""

    # Pairwise concordance
    html += """
    <div class="container">
        <h2>Pairwise Method Concordance</h2>
        <p>Agreement between pairs of methods (2-field resolution)</p>
        <table>
            <thead>
                <tr>
                    <th>Method Pair</th>
                    <th>Comparisons Made</th>
                    <th>Agreements (2-field)</th>
                    <th>Concordance Rate</th>
                    <th>1-field Concordance</th>
                </tr>
            </thead>
            <tbody>
"""

    method_names = {
        'optitype': 'OptiType',
        'hlala': 'HLA-LA',
        'hlahd': 'HLA-HD',
        'hisat': 'HISAT-genotype',
        'dragen': 'Dragen 4.4'
    }

    for pair_name, pair_stats in sorted(stats['method_pairs'].items()):
        if pair_stats['compared'] > 0:
            concordance_2f = (pair_stats['agree_2field'] / pair_stats['compared']) * 100
            concordance_1f = (pair_stats['agree_1field'] / pair_stats['compared']) * 100

            m1, m2 = pair_name.split('_vs_')
            pair_display = f"{method_names[m1]} vs {method_names[m2]}"

            css_class = ''
            if concordance_2f >= 80:
                css_class = 'high-concordance'
            elif concordance_2f >= 60:
                css_class = 'medium-concordance'
            else:
                css_class = 'low-concordance'

            html += f"""                <tr class="{css_class}">
                    <td>{pair_display}</td>
                    <td>{pair_stats['compared']}</td>
                    <td>{pair_stats['agree_2field']}</td>
                    <td><strong>{concordance_2f:.1f}%</strong></td>
                    <td>{concordance_1f:.1f}%</td>
                </tr>
"""

    html += """            </tbody>
        </table>
    </div>
"""

    # By gene concordance
    html += """
    <div class="container">
        <h2>Concordance by Gene</h2>
        <table>
            <thead>
                <tr>
                    <th>Gene</th>
                    <th>Samples</th>
                    <th>4 Methods Available</th>
                    <th>3 Methods Available</th>
                    <th>2 Methods Available</th>
                    <th>High Confidence Calls*</th>
                </tr>
            </thead>
            <tbody>
"""

    for gene in sorted(stats['by_gene'].keys()):
        gene_stats = stats['by_gene'][gene]
        total_samples = gene_stats['total_samples']
        high_conf = gene_stats['full_agreement']
        high_conf_pct = (high_conf / total_samples * 100) if total_samples > 0 else 0

        css_class = ''
        if high_conf_pct >= 80:
            css_class = 'high-concordance'
        elif high_conf_pct >= 50:
            css_class = 'medium-concordance'
        elif high_conf_pct > 0:
            css_class = 'low-concordance'

        html += f"""                <tr class="{css_class}">
                    <td><strong>{gene}</strong></td>
                    <td>{total_samples}</td>
                    <td>{gene_stats['all_4_methods']}</td>
                    <td>{gene_stats['any_3_methods']}</td>
                    <td>{gene_stats['any_2_methods']}</td>
                    <td>{high_conf} ({high_conf_pct:.0f}%)</td>
                </tr>
"""

    html += """            </tbody>
        </table>
        <p><small>*High Confidence = All available methods agree at 2-field resolution</small></p>
    </div>
"""

    # Interpretation guide
    html += """
    <div class="container">
        <h2>Interpretation Guide</h2>
        <h3>Reliability Scoring:</h3>
        <ul>
            <li><strong class="high-concordance" style="padding: 3px 8px;">All 4 methods agree:</strong> Highest confidence - allele is very likely correct</li>
            <li><strong class="high-concordance" style="padding: 3px 8px;">3 methods agree:</strong> High confidence - allele is likely correct</li>
            <li><strong class="medium-concordance" style="padding: 3px 8px;">2 methods agree:</strong> Moderate confidence - requires manual review</li>
            <li><strong class="low-concordance" style="padding: 3px 8px;">Only 1 method calls:</strong> Low confidence - likely false positive or rare variant</li>
        </ul>

        <h3>Resolution Levels:</h3>
        <ul>
            <li><strong>2-field (e.g., A*02:01):</strong> Protein-level resolution - used for clinical/transplant matching</li>
            <li><strong>1-field (e.g., A*02):</strong> Allele group - broader classification</li>
        </ul>

        <h3>Recommendations:</h3>
        <ul>
            <li>Use calls where ≥3 methods agree for high-confidence results</li>
            <li>Manually review calls where only 2 methods agree</li>
            <li>Consider single-method calls as low confidence unless from HLA-HD or HISAT-genotype for class II genes</li>
            <li>OptiType only types A, B, C - agreement limited to these genes</li>
        </ul>
    </div>
"""

    html += """
</body>
</html>
"""

    with open(output_file, 'w') as f:
        f.write(html)

def main():
    comparison_dir = Path('comparison_results')
    output_file = 'hla_concordance_analysis.html'

    if not comparison_dir.exists():
        print(f"Error: {comparison_dir} not found")
        return

    # Parse all HTML files
    html_files = list(comparison_dir.glob('*_comparison.html'))

    if not html_files:
        print(f"No comparison HTML files found in {comparison_dir}")
        return

    print(f"Analyzing concordance across {len(html_files)} samples...")

    all_samples = []
    for html_file in sorted(html_files):
        print(f"  Parsing {html_file.name}...")
        data = parse_html_comparison(html_file)
        all_samples.append(data)

    # Calculate concordance
    print("\nCalculating concordance statistics...")
    stats = calculate_concordance(all_samples)

    # Generate report
    print("Generating concordance report...")
    generate_concordance_html(stats, output_file)

    # Print summary
    total = stats['overall']['total_comparisons']
    if total > 0:
        all4_pct = (stats['overall']['all_4_agree_2field'] / total) * 100
        any3_pct = (stats['overall']['any_3_agree_2field'] / total) * 100

        print("\n" + "="*60)
        print("CONCORDANCE SUMMARY")
        print("="*60)
        print(f"\nTotal comparisons: {total}")
        print(f"All 4 methods agree: {stats['overall']['all_4_agree_2field']} ({all4_pct:.1f}%)")
        print(f"3+ methods agree: {stats['overall']['any_3_agree_2field']} ({any3_pct:.1f}%)")
        print(f"High confidence calls: {stats['overall']['all_4_agree_2field'] + stats['overall']['any_3_agree_2field']} ({all4_pct + any3_pct:.1f}%)")

        print("\nPairwise Concordance:")
        for pair_name, pair_stats in sorted(stats['method_pairs'].items()):
            if pair_stats['compared'] > 0:
                concordance = (pair_stats['agree_2field'] / pair_stats['compared']) * 100
                print(f"  {pair_name:25s}: {concordance:5.1f}% ({pair_stats['agree_2field']}/{pair_stats['compared']})")

    print(f"\n✓ Done! Open {output_file} to view the full concordance analysis")

if __name__ == '__main__':
    main()
