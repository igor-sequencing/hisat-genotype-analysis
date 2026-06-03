#!/bin/bash

# Auto-build HISAT2 graph indices after extract_vars completes

EXTRACT_VARS_LOG="extract_vars_full_run.log"
BUILD_LOG="hisat2_build.log"
INDICES_DIR="indicies"

echo "=== Monitoring extract_vars process ===" | tee -a "$BUILD_LOG"
echo "Started: $(date)" | tee -a "$BUILD_LOG"

# Wait for extract_vars to complete
while ps aux | grep -v grep | grep "hisatgenotype_extract_vars.py --base hla --locus-list" > /dev/null; do
    echo "$(date): extract_vars still running... ($(wc -l < $EXTRACT_VARS_LOG) lines in log)" | tee -a "$BUILD_LOG"
    sleep 30
done

echo "$(date): extract_vars completed!" | tee -a "$BUILD_LOG"
echo "" | tee -a "$BUILD_LOG"

# Check if the necessary files were created
echo "=== Checking for required files ===" | tee -a "$BUILD_LOG"
REQUIRED_FILES=(
    "$INDICES_DIR/hla_backbone.fa"
    "$INDICES_DIR/hla.index.snp"
    "$INDICES_DIR/hla.haplotype"
)

ALL_EXIST=true
for file in "${REQUIRED_FILES[@]}"; do
    if [ -f "$file" ]; then
        echo "✓ Found: $file ($(stat -c%s "$file" | numfmt --to=iec-i --suffix=B))" | tee -a "$BUILD_LOG"
    else
        echo "✗ Missing: $file" | tee -a "$BUILD_LOG"
        ALL_EXIST=false
    fi
done

if [ "$ALL_EXIST" = false ]; then
    echo "" | tee -a "$BUILD_LOG"
    echo "ERROR: Required files missing. Cannot build graph indices." | tee -a "$BUILD_LOG"
    echo "Check extract_vars_full_run.log for errors." | tee -a "$BUILD_LOG"
    exit 1
fi

# Count how many genes are in the backbone
GENE_COUNT=$(grep -c "^>" "$INDICES_DIR/hla_backbone.fa")
echo "" | tee -a "$BUILD_LOG"
echo "Found $GENE_COUNT genes in hla_backbone.fa" | tee -a "$BUILD_LOG"

# Show the gene list
echo "Genes:" | tee -a "$BUILD_LOG"
grep "^>" "$INDICES_DIR/hla_backbone.fa" | sed 's/^>/  - /' | tee -a "$BUILD_LOG"

echo "" | tee -a "$BUILD_LOG"
echo "=== Building HISAT2 graph indices ===" | tee -a "$BUILD_LOG"
echo "Started: $(date)" | tee -a "$BUILD_LOG"

# Run hisat2-build
hisat2-build \
    -p 4 \
    --snp "$INDICES_DIR/hla.index.snp" \
    --haplotype "$INDICES_DIR/hla.haplotype" \
    "$INDICES_DIR/hla_backbone.fa" \
    "$INDICES_DIR/hla.graph" 2>&1 | tee -a "$BUILD_LOG"

BUILD_EXIT_CODE=${PIPESTATUS[0]}

echo "" | tee -a "$BUILD_LOG"
echo "Completed: $(date)" | tee -a "$BUILD_LOG"
echo "Exit code: $BUILD_EXIT_CODE" | tee -a "$BUILD_LOG"

# Check if graph files were created
echo "" | tee -a "$BUILD_LOG"
echo "=== Verifying graph index files ===" | tee -a "$BUILD_LOG"

GRAPH_FILES_CREATED=0
for i in {1..8}; do
    GRAPH_FILE="$INDICES_DIR/hla.graph.$i.ht2"
    if [ -f "$GRAPH_FILE" ]; then
        SIZE=$(stat -c%s "$GRAPH_FILE" | numfmt --to=iec-i --suffix=B)
        echo "✓ Created: hla.graph.$i.ht2 ($SIZE)" | tee -a "$BUILD_LOG"
        ((GRAPH_FILES_CREATED++))
    else
        echo "✗ Missing: hla.graph.$i.ht2" | tee -a "$BUILD_LOG"
    fi
done

echo "" | tee -a "$BUILD_LOG"
if [ $GRAPH_FILES_CREATED -eq 8 ]; then
    echo "SUCCESS! All 8 graph index files created." | tee -a "$BUILD_LOG"
    echo "You can now run HLA genotyping with all $GENE_COUNT genes." | tee -a "$BUILD_LOG"
else
    echo "ERROR: Only $GRAPH_FILES_CREATED/8 graph files created." | tee -a "$BUILD_LOG"
    echo "Check hisat2_build.log for errors." | tee -a "$BUILD_LOG"
    exit 1
fi

echo "" | tee -a "$BUILD_LOG"
echo "=== Summary ===" | tee -a "$BUILD_LOG"
echo "Total indices size: $(du -sh $INDICES_DIR | cut -f1)" | tee -a "$BUILD_LOG"
echo "Log file: $BUILD_LOG" | tee -a "$BUILD_LOG"
