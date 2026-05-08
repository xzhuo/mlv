# Configuration file for Kraken2 installation and database setup. 
# Use the full path for KRAKEN2_INSTALL_DIR and KRAKEN_DB to avoid issues with relative paths. 
# THREADS is only used for the kraken2_run.sh.
# Make sure GENOME_DIR matches the one used in Snakefile.
KRAKEN2_INSTALL_DIR=~/kraken2
KRAKEN_DB=~/kraken2_db
THREADS=4
GENOME_DIR=genome

SRR=( \
	DRR349010 \
	# DRR349011 \
	# DRR349012 \
 	# DRR349013 \
 	# DRR349014 \
 	# DRR349015 \
	);