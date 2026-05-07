#!/bin/bash

source config.sh

# create the combined genome:
mkdir -p $GENOME_DIR
cd $GENOME_DIR
ln -s $KRAKEN_DB/library/9606/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna .
ln -s $KRAKEN_DB/library/10090/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.fna .
cat GCF_009914755.1_T2T-CHM13v2.0_genomic.fna GCF_000001635.27_GRCm39_genomic.fna > combined_genome.fna
bwa index combined_genome.fna

# index and create chromosome size files:
samtools faidx GCF_009914755.1_T2T-CHM13v2.0_genomic.fna
samtools faidx GCF_000001635.27_GRCm39_genomic.fna
perl -lne '@F=split /\||\s+/;print join("\t",@F[2,3])' GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.fai > GCF_009914755.1_T2T-CHM13v2.0_genomic.chrom_size
perl -lne '@F=split /\||\s+/;print join("\t",@F[2,3])' GCF_000001635.27_GRCm39_genomic.fna.fai > GCF_000001635.27_GRCm39_genomic.chrom_size
cd ../
