#!/bin/bash

source config.sh


# conda init;
# conda create --file envs/kmer.yml;
# conda activate kmer;

mkdir -p fastq;
cd fastq;
for id in ${SRR[@]};do
	mkdir -p ${id};
	aws s3 sync s3://sra-pub-run-odp/sra/${id} ${id} --no-sign-request;
	if [ -f ${id}/${id} ]; then
		fasterq-dump ${id}/${id} --threads $THREADS;
	else
		fasterq-dump ${id} --threads $THREADS;
	fi
	if [ -f ${id}_1.fastq ]; then
		rm -rf ${id};
		$KRAKEN2_INSTALL_DIR/kraken2 --db $KRAKEN_DB --threads $THREADS --output ${id}_1.txt ${id}_1.fastq;
		$KRAKEN2_INSTALL_DIR/kraken2 --db $KRAKEN_DB --threads $THREADS --output ${id}_2.txt ${id}_2.fastq;
		paste <(cut -f2,3,4,5 ${id}_1.txt) <(cut -f3,5 ${id}_2.txt) | perl -lne '@F=split /\t/;next unless ($F[1] == 10090 || $F[1] == 9606) && ($F[4] == 10090 || $F[4] ==9606); print $_ if $F[1] ne $F[4]' > ${id}_mosaic.txt;
		parallel "seqtk subseq ${id}_{}.fastq <(cut -f1 ${id}_mosaic.txt) > ${id}_mosaic_{}.fastq" ::: 1 2;
		python3 ../scripts/fqDustmasker.py -1 ${id}_mosaic_1.fastq -2 ${id}_mosaic_2.fastq -o ${id}_mosaic -t $THREADS;
		echo $(grep -c '^@' ${id}_mosaic_dust_1.fastq) $(wc -l ${id}_mosaic.txt) $(wc -l ${id}_1.txt) >> summary.txt;
		rm ${id}_1.fastq;
		rm ${id}_1.txt;
		# rm ${id}_mosaic_1.fastq;
		rm ${id}_2.fastq;
		rm ${id}_2.txt;
		# rm ${id}_mosaic_2.fastq;
	else
		echo "${id} dump fail!";
	fi
done;
