#!/usr/bin/env bash

# Step 1. Download the latest GENCODE files
# Visit: https://www.gencodegenes.org/human/
# Fasta files > Transcript sequences > Right Click on Fasta > Copy link address > Paste here 
TRANS_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.pc_transcripts.fa.gz"


# Visit: https://www.gencodegenes.org/human/
# Fasta files > Genome sequence, primary assembly (GRCh38) > Right Click on Fasta > Copy link address > Paste here 
GENOME_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/GRCh38.primary_assembly.genome.fa.gz"

set -euo pipefail
mkdir -p ../../../inputs/reference

wget -c "${TRANS_URL}" -O ../../../inputs/reference/gencode.v49.pc_transcripts.fa.gz
wget -c "${GENOME_URL}" -O ../../../inputs/reference/GRCh38.primary_assembly.genome.fa.gz

# Step 2. Create the decoy list
grep '^>' <(gunzip -c ../../../inputs/reference/GRCh38.primary_assembly.genome.fa.gz) | 
cut -d ' ' -f 1 > ../../../inputs/reference/decoys.txt
sed -i -e 's/>//g' ../../../inputs/reference/decoys.txt

# Step 3. Combine transcriptome and genome FASTAs
cat ../../../inputs/reference/gencode.v49.pc_transcripts.fa.gz \
    ../../../inputs/reference/GRCh38.primary_assembly.genome.fa.gz \
    > ../../../inputs/reference/transcripts_and_decoys.fa.gz

# Step 4. Build the Salmon index
salmon index \
  -t ../../../inputs/reference/transcripts_and_decoys.fa.gz \
  -d ../../../inputs/reference/decoys.txt \
  -p 30 \
  -i ../../../inputs/reference/human_salmon_index \
  --gencode

echo "All files and outputs saved in the existing ../../../input/reference directory"
