FILES = glob_wildcards('fastq/{sra}_mosaic_dust_1.fastq')
SRAS = FILES.sra

# provided in the repo:
RLTR4_MM = "genome/RLTR4_Mm.mm39.ncbi.bed"

# please make sure the dir matches the GENOME_DIR in the config.sh file.
COMBINED = "genome/combined_genome.fna"
MOUSE_SIZE = "genome/GCF_000001635.27_GRCm39_genomic.chrom_size"
HUMAN_SIZE = "genome/GCF_009914755.1_T2T-CHM13v2.0_genomic.chrom_size"

rule all:
  input:
    expand("fastp_fq/{sra}_fastp_{n}.fastq", sra=SRAS, n=[1, 2]),
    expand("fastp_fq/{sra}_fastp.json", sra=SRAS),
    expand("fastp_fq/{sra}_fastp.html", sra=SRAS),
    expand("bam/{sra}_fastp.sort.bam", sra=SRAS),
    expand("bam/{sra}_picard.txt", sra=SRAS),
    expand("bam/{sra}_dedup.bam", sra=SRAS),
    "all_dedup.bam",
    expand("bedpe/{sra}_dedup.bedpe", sra=SRAS),
    expand("bed/{sra}.{ref}.bed", sra=SRAS, ref=["human", "mouse"]),
    expand("all_dedup.{ltr}.fisher.txt", ltr=["RLTR4_Mm"]),
    expand("all_dedup.{ltr}.list", ltr=["RLTR4_Mm"])

rule cat_list:
  input:
    expand("bed/{sra}.{ltr}.list", sra=SRAS,allow_missing=True)
  output:
    "all_dedup.{ltr}.list"
  conda:
    "envs/align.yml"
  shell:
    "cat {input} > {output}"

rule cat_fisher:
  input:
    expand("bed/{sra}.{ltr}.fisher.txt", sra=SRAS,allow_missing=True)
  output:
    "all_dedup.{ltr}.fisher.txt"
  conda:
    "envs/align.yml"
  shell:
    "cat {input} > {output}"

rule cat_bam:
  input:
    expand("bam/{sra}_dedup.bam", sra=SRAS)
  output:
    "all_dedup.bam"
  conda:
    "envs/align.yml"
  shell:
    "samtools cat -o {output} {input}"


rule fisher:
  input:
    bed = "bed/{sra}.mouse.bed",
    rmsk = RLTR4_MM,
    gsize = MOUSE_SIZE
  output:
    fisher = "bed/{sra}.RLTR4_Mm.fisher.txt",
    lst = "bed/{sra}.RLTR4_Mm.list"
  params:
    cmd = r"""{$num = $F[-1] if /Number of overlaps/;($left,$right,$two,$ratio) = @F if EOF}END{print join("\t",$s,$te,$num,$ratio,$two,$left,$right)}"""
  conda:
    "envs/align.yml"
  shell:
    """
    bedtools intersect -a {input.bed} -b {input.rmsk} -u |cut -f4 > {output.lst};
    bedtools fisher -a {input.bed} -b {input.rmsk} -g {input.gsize} | perl -slane {params.cmd:q} -- -s={wildcards.sra} -te="RLTR4_Mm" > {output.fisher}
    """

rule split_bedpe: # splite bedped to a human bed and a mouse bed file. use the raw string so I don't have to escape everything anymore!
  input:
    "bedpe/{sra}_dedup.bedpe"
  params:
    cmd = r"""{($kraken,$ref1,$chr1)=split /\|/, $F[0];
         ($kraken,$ref2,$chr2)=split /\|/,$F[3];
         push @{$hash{$ref1}}, [ $chr1,$F[1],$F[2],$F[6],$F[7],$F[8] ];
         push @{$hash{$ref2}},[ $chr2,$F[4],$F[5],$F[6],$F[7],$F[9] ]
         }END{
           @mouse = sort {$a->[0] cmp $b->[0] || $a->[1] <=> $b->[1]} @{$hash{"10090"}};
           @human = sort {$a->[0] cmp $b->[0] || $a->[1] <=> $b->[1]} @{$hash{"9606"}};
           open(A, ">", "bed/$sra.mouse.bed");
           foreach $i (@mouse){print A join("\t",@{$i})};
           close A;
           open(B, ">", "bed/$sra.human.bed");
           foreach $i (@human){print B join("\t",@{$i})};
           close B;}"""

  output:
    human = "bed/{sra}.human.bed",
    mouse = "bed/{sra}.mouse.bed"
  conda:
    "envs/align.yml"
  shell:
    "perl -slane {params.cmd:q} -- -sra={wildcards.sra} {input}"

rule bedpe:
  input:
    "bam/{sra}_dedup.bam"
  output:
    "bedpe/{sra}_dedup.bedpe"
  params:
    cmd = r"""next if $F[0] eq "." || $F[3] eq ".";($kraken,$ref1,$chr1)=split /\|/, $F[0];($kraken,$ref2,$chr2)=split /\|/,$F[3];next if $ref1 eq $ref2; print $_"""
  conda:
    "envs/align.yml"
  shell:
    "samtools view -b -F 0X800 -F 0X400 {input} |samtools collate -O - | bedtools bamtobed -bedpe -i stdin| perl -lane {params.cmd:q} > {output}"

rule markdup:
  input:
    "bam/{sra}_fastp.sort.bam"
  output:
    bam = "bam/{sra}_dedup.bam",
    m = "bam/{sra}_picard.txt"
  conda:
    "envs/align.yml"
  shell:
    "picard MarkDuplicates -I {input} -O {output.bam} -M {output.m}"

rule bwa:
  input:
    expand("fastp_fq/{sra}_fastp_{n}.fastq", n=[1, 2], allow_missing=True)
  output:
    "bam/{sra}_fastp.sort.bam"
  threads:
    workflow.cores
  conda:
    "envs/align.yml"
  shell:
    """bwa mem -5SP -t {threads} {COMBINED} {input[0]} {input[1]} |samtools addreplacerg -r "@RG\tID:{wildcards.sra}\tSM:{wildcards.sra}\tPL:Illumina\tLB:{wildcards.sra}" - | samtools sort -o {output}"""

rule fastp:
  input:
    expand("fastq/{sra}_mosaic_dust_{n}.fastq", n=[1, 2], allow_missing=True)
  output:
    fq = expand("fastp_fq/{sra}_fastp_{n}.fastq", n=[1, 2], allow_missing=True),
    json = "fastp_fq/{sra}_fastp.json",
    html = "fastp_fq/{sra}_fastp.html"
  conda:
    "envs/align.yml"
  shell:
    """fastp --dedup -i {input[0]} -I {input[1]} -j {output.json} -h {output.html} -o {output.fq[0]} -O {output.fq[1]}"""
