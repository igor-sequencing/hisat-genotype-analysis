#!/usr/bin/env python3
"""
Merge HLA typing results from multiple tools.
Uses organized input structure from hla_typing_inputs/ directory.
Create separate HTML file for each sample.
Rows: genes, Columns: methods
"""

import os
import re
import csv
from pathlib import Path
from collections import defaultdict

def parse_hisat_report(report_path):
    """Parse HISAT-genotype report and extract top ranked alleles per gene."""
    results = {}

    with open(report_path, 'r') as f:
        content = f.read()

    # Find all ranked sections - get top 2 ranked alleles
    pattern = r'(\d+) ranked ([A-Z0-9]+)\*([^\s]+) \(abundance: ([\d.]+)%\)'
    matches = re.findall(pattern, content)

    for rank, gene, allele, abundance in matches:
        if gene not in results:
            results[gene] = []

        # Only keep top 2 ranked
        if int(rank) <= 2:
            results[gene].append({
                'rank': int(rank),
                'allele': f"HLA-{gene}*{allele}",
                'abundance': float(abundance)
            })

    return results

def parse_hlahd_result(result_path):
    """Parse HLA-HD final result file."""
    results = {}

    with open(result_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            parts = line.split('\t')
            if len(parts) >= 3:
                gene = parts[0].strip()
                allele1 = parts[1].strip()
                allele2 = parts[2].strip()

                results[gene] = []
                if allele1 and allele1 != '-' and allele1 != 'Not typed':
                    results[gene].append(allele1)
                if allele2 and allele2 != '-' and allele2 != 'Not typed':
                    results[gene].append(allele2)

    return results

def parse_optitype_result(result_path):
    """Parse OptiType TSV result file."""
    results = {}

    if not Path(result_path).exists():
        return results

    # Read the TSV file
    with open(result_path, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            # OptiType provides A1, A2, B1, B2, C1, C2
            for gene in ['A', 'B', 'C']:
                alleles = []
                for i in ['1', '2']:
                    col = f'{gene}{i}'
                    if col in row and row[col]:
                        allele = row[col].strip()
                        if allele:
                            # Format as HLA-A*02:01
                            alleles.append(f'HLA-{allele}')

                if alleles:
                    results[gene] = alleles

    return results

def parse_hlala_result(bestguess_path):
    """Parse HLA-LA bestguess file."""
    results = {}

    if not Path(bestguess_path).exists():
        return results

    # Read the file
    with open(bestguess_path, 'r') as f:
        lines = f.readlines()

    # Skip header line
    for line in lines[1:]:
        line = line.strip()
        if not line:
            continue

        parts = line.split('\t')
        if len(parts) >= 3:
            gene = parts[0].strip()
            allele = parts[2].strip()

            # Remove G suffix and N suffix for cleaner display
            allele = allele.replace('G', '').strip()

            if gene not in results:
                results[gene] = []

            # Format as HLA-A*02:01:01
            if not allele.startswith('HLA-'):
                allele = f'HLA-{allele}'

            results[gene].append(allele)

    return results

def parse_dragen_result(tsv_path):
    """Parse Dragen HLA TSV result file."""
    results = {}

    if not Path(tsv_path).exists():
        return results

    # Read the TSV file
    with open(tsv_path, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            gene = row.get('gene', '').strip()
            if not gene or gene in ['HFE', 'Y', 'R']:  # Skip non-HLA genes
                continue

            num_alleles = row.get('num_alleles', '0').strip()
            if num_alleles == '0':
                continue

            alleles = []
            # Get allele 1
            allele1 = row.get('allele_1', '').strip()
            if allele1 and allele1 != 'NA':
                # Format as HLA-A*02:01
                if not allele1.startswith('HLA-'):
                    allele1 = f'HLA-{allele1}'
                alleles.append(allele1)

            # Get allele 2 if present
            allele2 = row.get('allele_2', '').strip()
            if allele2 and allele2 != 'NA':
                if not allele2.startswith('HLA-'):
                    allele2 = f'HLA-{allele2}'
                alleles.append(allele2)

            if alleles:
                results[gene] = alleles

    return results

def collect_sample_data(sample_name, input_base_dir):
    """Collect data from all tools for a sample."""
    data = {
        'sample': sample_name,
        'hisat': {},
        'hlahd': {},
        'optitype': {},
        'hlala': {},
        'dragen': {}
    }

    base_path = Path(input_base_dir)

    # Get HISAT data
    hisat_file = base_path / 'hisat' / f'{sample_name}.report'
    if hisat_file.exists():
        data['hisat'] = parse_hisat_report(hisat_file)

    # Get HLA-HD data
    hlahd_file = base_path / 'hlahd' / f'{sample_name}_final.result.txt'
    if hlahd_file.exists():
        data['hlahd'] = parse_hlahd_result(hlahd_file)

    # Get OptiType data
    optitype_file = base_path / 'optitype' / f'{sample_name}_result.tsv'
    if optitype_file.exists():
        data['optitype'] = parse_optitype_result(optitype_file)

    # Get HLA-LA data
    hlala_file = base_path / 'hlala' / f'{sample_name}_bestguess.txt'
    if hlala_file.exists():
        data['hlala'] = parse_hlala_result(hlala_file)

    # Get Dragen data
    dragen_file = base_path / 'dragen' / f'{sample_name}_HLA.tsv'
    if dragen_file.exists():
        data['dragen'] = parse_dragen_result(dragen_file)

    return data

def generate_sample_html(data, output_file):
    """Generate HTML table for a single sample comparing methods."""
    sample_name = data['sample']

    # Get all genes from all sources
    all_genes = set()
    all_genes.update(data['hisat'].keys())
    all_genes.update(data['hlahd'].keys())
    all_genes.update(data['optitype'].keys())
    all_genes.update(data['hlala'].keys())
    all_genes.update(data['dragen'].keys())
    all_genes = sorted(all_genes)

    html = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>HLA Typing Results - {sample_name}</title>
    <style>
        body {{
            font-family: Arial, sans-serif;
            margin: 20px;
            background-color: #f5f5f5;
        }}
        h1 {{
            color: #333;
            text-align: center;
        }}
        .sample-info {{
            text-align: center;
            margin-bottom: 20px;
            font-size: 18px;
            color: #666;
        }}
        .table-container {{
            overflow-x: auto;
            background-color: white;
            padding: 20px;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            max-width: 1400px;
            margin: 0 auto;
        }}
        table {{
            border-collapse: collapse;
            width: 100%;
        }}
        th, td {{
            border: 1px solid #ddd;
            padding: 12px;
            text-align: left;
        }}
        th {{
            background-color: #4CAF50;
            color: white;
            font-weight: bold;
        }}
        th.gene-col {{
            background-color: #2196F3;
            width: 100px;
        }}
        th.method-col {{
            background-color: #FF9800;
            width: 280px;
        }}
        tr:nth-child(even) {{
            background-color: #f9f9f9;
        }}
        tr:hover {{
            background-color: #f1f1f1;
        }}
        .allele {{
            display: block;
            margin: 3px 0;
            padding: 2px 5px;
            background-color: #e3f2fd;
            border-radius: 3px;
            font-family: monospace;
            font-size: 11px;
        }}
        .optitype-allele {{
            background-color: #f3e5f5;
        }}
        .hlala-allele {{
            background-color: #e8f5e9;
        }}
        .hlahd-allele {{
            background-color: #fff3e0;
        }}
        .dragen-allele {{
            background-color: #fce4ec;
        }}
        .rank-info {{
            color: #666;
            font-size: 10px;
        }}
        .no-data {{
            color: #999;
            font-style: italic;
        }}
    </style>
</head>
<body>
    <h1>HLA Typing Results Comparison</h1>
    <div class="sample-info">Sample: <strong>{sample_name}</strong></div>
    <div class="table-container">
        <table>
            <thead>
                <tr>
                    <th class="gene-col">Gene</th>
                    <th class="method-col">OptiType</th>
                    <th class="method-col">HLA-LA</th>
                    <th class="method-col">HLA-HD</th>
                    <th class="method-col">HISAT-genotype</th>
                    <th class="method-col">Dragen 4.4</th>
                </tr>
            </thead>
            <tbody>
"""

    # Add rows for each gene
    for gene in all_genes:
        html += f"                <tr>\n"
        html += f"                    <td><strong>{gene}</strong></td>\n"

        # OptiType column - sort alphabetically
        html += "                    <td>"
        if gene in data['optitype'] and data['optitype'][gene]:
            for allele in sorted(data['optitype'][gene]):
                html += f'<span class="allele optitype-allele">{allele}</span>'
        else:
            html += '<span class="no-data">No data</span>'
        html += "</td>\n"

        # HLA-LA column - sort alphabetically
        html += "                    <td>"
        if gene in data['hlala'] and data['hlala'][gene]:
            for allele in sorted(data['hlala'][gene]):
                html += f'<span class="allele hlala-allele">{allele}</span>'
        else:
            html += '<span class="no-data">No data</span>'
        html += "</td>\n"

        # HLA-HD column - sort alphabetically
        html += "                    <td>"
        if gene in data['hlahd'] and data['hlahd'][gene]:
            for allele in sorted(data['hlahd'][gene]):
                html += f'<span class="allele hlahd-allele">{allele}</span>'
        else:
            html += '<span class="no-data">No data</span>'
        html += "</td>\n"

        # HISAT-genotype column - sort alphabetically, preserve rank as attribute
        html += "                    <td>"
        if gene in data['hisat'] and data['hisat'][gene]:
            for allele_info in sorted(data['hisat'][gene], key=lambda x: x['allele']):
                allele = allele_info['allele']
                abundance = allele_info['abundance']
                rank = allele_info['rank']
                html += f'<span class="allele" data-rank="{rank}">{allele} <span class="rank-info">(rank: {rank}, {abundance}%)</span></span>'
        else:
            html += '<span class="no-data">No data</span>'
        html += "</td>\n"

        # Dragen column - sort alphabetically
        html += "                    <td>"
        if gene in data['dragen'] and data['dragen'][gene]:
            for allele in sorted(data['dragen'][gene]):
                html += f'<span class="allele dragen-allele">{allele}</span>'
        else:
            html += '<span class="no-data">No data</span>'
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

    print(f"Generated: {output_file}")

def main():
    input_base_dir = 'hla_typing_inputs'
    output_dir = 'comparison_results'

    # Create output directory
    Path(output_dir).mkdir(exist_ok=True)

    # Get all samples from hisat directory
    hisat_path = Path(input_base_dir) / 'hisat'
    if not hisat_path.exists():
        print(f"Error: {hisat_path} not found")
        return

    samples = []
    for report_file in sorted(hisat_path.glob('*.report')):
        sample_name = report_file.stem
        samples.append(sample_name)

    print(f"Processing {len(samples)} samples from {input_base_dir}...")

    # Process each sample
    for sample_name in samples:
        print(f"\nProcessing {sample_name}...")

        # Collect data for this sample
        data = collect_sample_data(sample_name, input_base_dir)

        # Generate HTML file
        output_file = Path(output_dir) / f'{sample_name}_comparison.html'
        generate_sample_html(data, output_file)

    print(f"\n✓ Done! Generated {len(samples)} comparison files in '{output_dir}/' directory")
    print(f"  Open any HTML file to view the comparison")

if __name__ == '__main__':
    main()
