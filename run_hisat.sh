# Use system samtools (upgraded to 1.22)
samtools view -b -M $1 6:28500000-33500000 > $2/mhc.bam
samtools sort -n -o $2/mhc.namesort.bam $2/mhc.bam

samtools fastq \
    -1 $2/mhc_R1.fastq \
    -2 $2/mhc_R2.fastq \
    -0 /dev/null -s /dev/null \
    $2/mhc.namesort.bam

cd /mnt/data/refData/hla-tools/hisat-genotype
# Set PYTHONPATH so hisatgenotype can find its modules
export PYTHONPATH=/mnt/data/refData/hla-tools/hisat-genotype/hisatgenotype_modules:$PYTHONPATH
# Add hisat2 to PATH
export PATH=/mnt/data/refData/hla-tools/hisat-genotype/hisat2:$PATH
# Reduced --pp from 20 to 8 to prevent hisat2 crashes
./hisatgenotype \
  -z indicies \
  --base hla \
  --pp 8 \
  -1 $2/mhc_R1.fastq \
  -2 $2/mhc_R2.fastq \
  --out-dir $2
