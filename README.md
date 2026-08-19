# Kmeric: kmer-based chimeric reads identifier

Murine Leukemia virus (MLV) may infect human cells. To efficiently identify MLV integration within human genome, we developed the kmeric pipeline to quickly identify human-MLV chimeric reads from patient-derived xenograft (PDX) paired-end sequencing libraries using Kraken2.

The pipeline can be modified to investigate other viral integration and cross-species transmission events.

## Pipeline overview

* Build a custom (human-mouse) Kraken2 datase.
* Download fastq files from SRA.
* Run Kraken2 on either end of paired fastq files separately.
* Retain only human-mouse chimeric paired-end reads.
* Remove reads with low-complexity or low-quality with dustMasker and fastp.
* Align remaining high quality human-mouse chimeric reads to a combined human-mouse reference genome.
* Catalog and calculate the enrichment of human-MLV chimeric reads among all mapped human-mouse chimeric reads.

## Manual

The pipeline can be installed on Linux x64 platform with conda.

1. Download the git repository:

```sh
git clone https://github.com/xzhuo/kmeric.git
```

Before started, you can modify the `config.sh` file to setup which SRA accessions to be used for the screening. You can also change Kraken and genome directory.

2. Install kraken and setup human, mouse kmer database:

```sh
bash kraken2_install.sh
```

3. Created two env (`kmer` and `align`) from yml files:

```sh
conda env create --file envs/kmer.yml
conda env create --file envs/align.yml
```

4. Run kraken2 with conda environment `kmer`:

```sh
conda run -n kmer bash kraken2_run.sh
```

5. Prepare necessary genome files with conda environment `align`:

```sh
conda run -n align bash genome.sh
```

6. Run the Snakemake pipeline (alignment and enrichment analysis) with conda environment `align`:

```sh
conda run -n align snakemake -j 4
```

## Output files

* all_dedup.bam
all deduplicated human-mouse chimeric reads aligned to the CHM13-mm39 genome in bam format.

* all_dedup.RLTR4_Mm.fisher.txt
The enrichment of human-MLV chimeric reads among all human-mouse chimeric reads calculated by Bedtools fisher.

    Columns in the file:

    > Column1: accession.
    >
    > Column2: should be "RLTR4_Mm" when testing the enrichment of MLV.
    >
    > Column3: number of human-MLV chimeric reads in each paired fastq files.
    >
    > Column4: enrichment ratio.
    >
    > Column5: two-tail P-value.
    >
    > Column6: left-tail P-value.
    >
    > Column7: right-tail P-value.

* all_dedup.RLTR4_Mm.list:
A list of identified human-MLV chimeric reads
