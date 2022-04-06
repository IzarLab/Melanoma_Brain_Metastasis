#!/bin/bash

### RNA velocity: Generate loom file (velocyto)

pat=$1

### Sync in from AWS
# for sc
if [[ "$pat" == *"_sc"* ]]; then
  aws s3 sync s3://scrna-seq/cellranger/v6.1.1/${pat}_CD45pos/ data/${pat}/ --exclude '*' \
   --include 'possorted_genome_bam*' --quiet
fi

# for sn
if [[ "$pat" == *"_sn"* ]]; then
  aws s3 sync s3://snrna-seq/cellranger/v6.1.1/${pat}/ data/${pat}/ --exclude '*' \
   --include 'possorted_genome_bam*' --quiet
fi

### Run velocyto
velocyto run data/${pat}/possorted_genome_bam.bam \
data/RNA_velocity/Homo_sapiens.GRCh38.102.gtf \
-m data/RNA_velocity/hm10_rmsk.gtf \
-@ 5 \
--sampleid ${pat} \
--outputfolder data/${pat} \
-vvv
