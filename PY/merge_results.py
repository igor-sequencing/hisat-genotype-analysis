#!/usr/bin/env python3
"""
Merge HISAT-genotype results with estimation results.
Create separate HTML file for each kit/sample.
Rows: genes, Columns: methods (estimation vs hisat-genotype)
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

def parse_estimation_result(result_path):
    """Parse estimation final result file."""
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

def parse_optitype_result(optitype_dir, sample_name):
    """Parse OptiType TSV result file."""
    results = {}

    # Find the result TSV file
    sample_path = Path(optitype_dir) / sample_name
    if not sample_path.exists():
        return results

    # Find the most recent result file
    tsv_files = list(sample_path.rglob('*_result.tsv'))
    if not tsv_files:
        return results

    # Read the TSV file
    with open(tsv_files[0], 'r') as f:
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

def parse_hlala_result(hlala_dir, sample_name):
    """Parse HLA-LA bestguess file."""
    results = {}

    # Find the bestguess file
    bestguess_file = Path(hlala_dir) / sample_name / 'hla' / 'R1_bestguess_G.txt'
    if not bestguess_file.exists():
        return results

    # Read the file
    with open(bestguess_file, 'r') as f:
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

def collect_sample_data(sample_name, hisat_dir, estimation_dir, optitype_dir, hlala_dir):
    """Collect HISAT, estimation, OptiType, and HLA-LA data for a sample."""
    data = {
        'sample': sample_name,
        'hisat': {},
        'estimation': {},
        'optitype': {},
        'hlala': {}
    }

    # Get HISAT data
    hisat_path = Path(hisat_dir) / sample_name
    if hisat_path.exists():
        report_files = list(hisat_path.glob('*.report'))
        if report_files:
            data['hisat'] = parse_hisat_report(report_files[0])

    # Get estimation data
    estimation_path = Path(estimation_dir) / sample_name / 'result' / f'{sample_name}_final.result.txt'
    if estimation_path.exists():
        data['estimation'] = parse_estimation_result(estimation_path)

    # Get OptiType data
    data['optitype'] = parse_optitype_result(optitype_dir, sample_name)

    # Get HLA-LA data
    data['hlala'] = parse_hlala_result(hlala_dir, sample_name)

    return data

def generate_sample_html(data, output_file):
    """Generate HTML table for a single sample comparing methods."""
    sample_name = data['sample']

    # Get all genes from all sources
    all_genes = set()
    all_genes.update(data['hisat'].keys())
    all_genes.update(data['estimation'].keys())
    all_genes.update(data['optitype'].keys())
    all_genes.update(data['hlala'].keys())
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
            max-width: 1200px;
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
            width: 350px;
        }}
        .optitype-allele {{
            background-color: #f3e5f5;
        }}
        .hlala-allele {{
            background-color: #e8f5e9;
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
        }}
        .estimation-allele {{
            background-color: #fff3e0;
        }}
        .rank-info {{
            color: #666;
            font-size: 11px;
        }}
        .no-data {{
            color: #999;
            font-style: italic;
        }}
        .matching {{
            background-color: #c8e6c9;
            font-weight: bold;
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

        # Estimation column - sort alphabetically
        html += "                    <td>"
        if gene in data['estimation'] and data['estimation'][gene]:
            for allele in sorted(data['estimation'][gene]):
                html += f'<span class="allele estimation-allele">{allele}</span>'
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
    hisat_dir = 'hisat.output'
    estimation_dir = '/mnt/data/igor-shared/HLA/estimation'
    optitype_dir = 'optitype_output'
    hlala_dir = '/igor-shared/HLA-LA/working'
    output_dir = 'comparison_results'

    # Create output directory
    Path(output_dir).mkdir(exist_ok=True)

    # Get all samples from hisat output
    hisat_path = Path(hisat_dir)
    if not hisat_path.exists():
        print(f"Error: {hisat_dir} not found")
        return

    samples = []
    for sample_dir in sorted(hisat_path.iterdir()):
        if sample_dir.is_dir():
            samples.append(sample_dir.name)

    print(f"Processing {len(samples)} samples...")

    # Process each sample
    for sample_name in samples:
        print(f"\nProcessing {sample_name}...")

        # Collect data for this sample
        data = collect_sample_data(sample_name, hisat_dir, estimation_dir, optitype_dir, hlala_dir)

        # Generate HTML file
        output_file = Path(output_dir) / f'{sample_name}_comparison.html'
        generate_sample_html(data, output_file)

    print(f"\n✓ Done! Generated {len(samples)} comparison files in '{output_dir}/' directory")
    print(f"  Open any HTML file to view the comparison")

if __name__ == '__main__':
    main()
