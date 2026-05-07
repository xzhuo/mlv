#!/bin/bash

source config.sh

# create the combined genome:
mkdir -p $GENOME_DIR
cat $KRAKEN_DB/library/9606/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna $KRAKEN_DB/library/10090/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.fna > $GENOME_DIR/combined_genome.fna