source config.sh

mkdir -p $KRAKEN2_INSTALL_DIR
mkdir -p $KRAKEN_DB
git clone https://github.com/DerrickWood/kraken2.git

kraken2/install_kraken2.sh $KRAKEN2_INSTALL_DIR
# You may have to specify the compiler if needed, for example with gcc-15:
# CC=gcc-15 CXX=g++-15 kraken2/install_kraken2.sh $KRAKEN2_INSTALL_DIR

$KRAKEN2_INSTALL_DIR/k2 download-library --db $KRAKEN_DB --library 10090 # Download mouse library
$KRAKEN2_INSTALL_DIR/k2 download-library --db $KRAKEN_DB --library 9606 # Download human library
$KRAKEN2_INSTALL_DIR/k2 download-taxonomy --db $KRAKEN_DB --skip-maps
$KRAKEN2_INSTALL_DIR/kraken2-build --build --db $KRAKEN_DB --threads 4 #Build the database
rm -rf kraken2

# create the combined genome:
cat $KRAKEN_DB/library/9606/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna $KRAKEN_DB/library/10090/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.fna > $COMBINED_GENOME