FILES = glob_wildcards('fastq/{sra}_mosaic_dust_1.fastq')
SRAS = FILES.sra

PICARD = "/opt/apps/picard/2.8.1/picard.jar"
HG38_MM10 = "/scratch/genomes/hg38_mm10/10xATAC/refdata-cellranger-atac-GRCh38-and-mm10-1.2.0/fasta/genome.fa"
MM10_RMSK = "/scratch/xzhuo/genomes/mm10/mm10.rmsk.class.lite.bed"
MLV = "/scratch/xzhuo/xenograft_fq/mlv.mm10.bed"
RLTR4_MM = "/scratch/xzhuo/genomes/mm10/RLTR4_Mm.mm10.bed"
RLTR6_MM = "/scratch/xzhuo/genomes/mm10/mm10.RLTR6_Mm.lite.bed"
IAPEY = "/scratch/xzhuo/genomes/mm10/mm10.IAPEY_LTR.lite.bed"
LTR5_HS = "/scratch/xzhuo/genomes/hg38/reAnno_hg38.LTR5_Hs.bed"
MM10_SIZE = "/scratch/genomes/mm10/mm10_lite.size"
HG38_SIZE = "/scratch/genomes/hg38/hg38_lite.size"
FISHER_PY = "/scratch/xzhuo/github/mlv/fisher.py"

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
    expand("bed/{sra}.{ref}.bed", sra=SRAS, ref=["hg38", "mm10"]),
    expand("bed/{sra}.mlv.count.txt", sra=SRAS),
    "all_dedup.mlv.count.txt",
    expand("all_dedup.{ltr}.fisher.txt", ltr=["RLTR4_Mm", "RLTR6_Mm", "IAPEY", "LTR5_Hs"]),
    expand("all_dedup.{ltr}.list", ltr=["RLTR4_Mm", "RLTR6_Mm", "IAPEY", "LTR5_Hs"])

rule cat_list:
  input:
    expand("bed/{sra}.{ltr}.list", sra=SRAS,allow_missing=True)
  output:
    "all_dedup.{ltr}.list"
  shell:
    "cat {input} > {output}"

rule cat_fisher:
  input:
    expand("bed/{sra}.{ltr}.fisher.txt", sra=SRAS,allow_missing=True)
  output:
    "all_dedup.{ltr}.fisher.txt"
  shell:
    "cat {input} > {output}"

rule cat_bam:
  input:
    expand("bam/{sra}_dedup.bam", sra=SRAS)
  output:
    "all_dedup.bam"
  shell:
    "ml samtools; samtools cat -o {output} {input}"

rule LTR5_fisher:
  input:
    bed = "bed/{sra}.hg38.bed",
    rmsk = LTR5_HS,
    gsize = HG38_SIZE
  output:
    fisher = "bed/{sra}.LTR5_Hs.fisher.txt",
    lst = "bed/{sra}.LTR5_Hs.list"
  params:
    cmd = r"""{$num = $F[-1] if /Number of overlaps/;($left,$right,$two,$ratio) = @F if EOF}END{print join("\t",$s,$te,$num,$ratio,$two,$left,$right)}"""
  shell:
    """
    /scratch/devtools/xzhuo/miniconda3/envs/pybedtools/bin/bedtools intersect -a {input.bed} -b {input.rmsk} -u |cut -f4 > {output.lst};
    /scratch/devtools/xzhuo/miniconda3/envs/pybedtools/bin/bedtools fisher -a {input.bed} -b {input.rmsk} -g {input.gsize} | perl -slane {params.cmd:q} -- -s={wildcards.sra} -te="LTR5_Hs" > {output.fisher}
    """
    # "python3 {FISHER_PY} -m LTR5_Hs -i {input.bed} -r {input.rmsk} -g {input.gsize} -o {output.fisher} -l {output.lst}"

rule IAPEY_fisher:
  input:
    bed = "bed/{sra}.mm10.bed",
    rmsk = IAPEY,
    gsize = MM10_SIZE
  output:
    fisher = "bed/{sra}.IAPEY.fisher.txt",
    lst = "bed/{sra}.IAPEY.list"
  params:
    cmd = r"""{$num = $F[-1] if /Number of overlaps/;($left,$right,$two,$ratio) = @F if EOF}END{print join("\t",$s,$te,$num,$ratio,$two,$left,$right)}"""
  shell:
    """
    /scratch/devtools/xzhuo/miniconda3/envs/pybedtools/bin/bedtools intersect -a {input.bed} -b {input.rmsk} -u |cut -f4 > {output.lst};
    /scratch/devtools/xzhuo/miniconda3/envs/pybedtools/bin/bedtools fisher -a {input.bed} -b {input.rmsk} -g {input.gsize} | perl -slane {params.cmd:q} -- -s={wildcards.sra} -te="IAPEY" > {output.fisher}
    """
    # "python3 {FISHER_PY} -m IAPEY -i {input.bed} -r {input.rmsk} -g {input.gsize} -o {output.fisher} -l {output.lst}"

rule RLTR6_fisher:
  input:
    bed = "bed/{sra}.mm10.bed",
    rmsk = RLTR6_MM,
    gsize = MM10_SIZE
  output:
    fisher = "bed/{sra}.RLTR6_Mm.fisher.txt",
    lst = "bed/{sra}.RLTR6_Mm.list"
  params:
    cmd = r"""{$num = $F[-1] if /Number of overlaps/;($left,$right,$two,$ratio) = @F if EOF}END{print join("\t",$s,$te,$num,$ratio,$two,$left,$right)}"""
  shell:
    """
    /scratch/devtools/xzhuo/miniconda3/envs/pybedtools/bin/bedtools intersect -a {input.bed} -b {input.rmsk} -u |cut -f4 > {output.lst};
    /scratch/devtools/xzhuo/miniconda3/envs/pybedtools/bin/bedtools fisher -a {input.bed} -b {input.rmsk} -g {input.gsize} | perl -slane {params.cmd:q} -- -s={wildcards.sra} -te="RLTR6_Mm" > {output.fisher}
    """
    # "python3 {FISHER_PY} -m RLTR6_Mm -i {input.bed} -r {input.rmsk} -g {input.gsize} -o {output.fisher} -l {output.lst}"

rule fisher:
  input:
    bed = "bed/{sra}.mm10.bed",
    rmsk = RLTR4_MM,
    gsize = MM10_SIZE
  output:
    fisher = "bed/{sra}.RLTR4_Mm.fisher.txt",
    lst = "bed/{sra}.RLTR4_Mm.list"
  params:
    cmd = r"""{$num = $F[-1] if /Number of overlaps/;($left,$right,$two,$ratio) = @F if EOF}END{print join("\t",$s,$te,$num,$ratio,$two,$left,$right)}"""
  # conda:
  #   "pybedtools"
  shell:
    """
    /scratch/devtools/xzhuo/miniconda3/envs/pybedtools/bin/bedtools intersect -a {input.bed} -b {input.rmsk} -u |cut -f4 > {output.lst};
    /scratch/devtools/xzhuo/miniconda3/envs/pybedtools/bin/bedtools fisher -a {input.bed} -b {input.rmsk} -g {input.gsize} | perl -slane {params.cmd:q} -- -s={wildcards.sra} -te="RLTR4_Mm" > {output.fisher}
    """
    # "python3 {FISHER_PY} -m RLTR4_Mm -i {input.bed} -r {input.rmsk} -g {input.gsize} -o {output.fisher} -l {output.lst}"

rule mlv_cat:
  input:
    expand("bed/{sra}.mlv.count.txt", sra=SRAS)
  output:
    "all_dedup.mlv.count.txt"
  conda:
    "env/align.yml"
  shell:
    "cat {input} > {output}"

rule mlv_count:
  input:
    bed = "bed/{sra}.mm10.bed",
    mlv = MLV
  params:
    cmd = r"""print join("\t",$sra,@F[3,6])"""
  output:
    "bed/{sra}.mlv.count.txt"
  conda:
    "env/align.yml"
  shell:
    "bedtools intersect -a {input.mlv} -b {input.bed} -c | perl -slane {params.cmd:q} -- -sra={wildcards.sra} > {output}"

rule split_bedpe: # splite bedped to a hg38 bed and a mm10 bed file. use the raw string so I don't have to escape everything anymore!
  input:
    "bedpe/{sra}_dedup.bedpe"
  params:
    cmd = r"""{($ref1,$chr1)=split /_/, $F[0], 2;
         ($ref2,$chr2)=split /_/,$F[3], 2;
         push @{$hash{$ref1}}, [ $chr1,$F[1],$F[2],$F[6],$F[7],$F[8] ];
         push @{$hash{$ref2}},[ $chr2,$F[4],$F[5],$F[6],$F[7],$F[9] ]
         }END{
           @mm10 = sort {$a->[0] cmp $b->[0] || $a->[1] <=> $b->[1]} @{$hash{"mm10"}};
           @hg38 = sort {$a->[0] cmp $b->[0] || $a->[1] <=> $b->[1]} @{$hash{"GRCh38"}};
           open(A, ">", "bed/$sra.mm10.bed");
           foreach $i (@mm10){print A join("\t",@{$i})};
           close A;
           open(B, ">", "bed/$sra.hg38.bed");
           foreach $i (@hg38){print B join("\t",@{$i})};
           close B;}"""

  output:
    hg38 = "bed/{sra}.hg38.bed",
    mm10 = "bed/{sra}.mm10.bed"
  conda:
    "env/align.yml"
  shell:
    "perl -slane {params.cmd:q} -- -sra={wildcards.sra} {input}"

rule bedpe:
  input:
    "bam/{sra}_dedup.bam"
  output:
    "bedpe/{sra}_dedup.bedpe"
  params:
    cmd = r"""next if $F[0] eq "." || $F[3] eq ".";($ref1,$chr1)=split /_/, $F[0], 2;($ref2,$chr2)=split /_/,$F[3], 2;next if $ref1 eq $ref2; next if length($chr1)>=6 || length($chr2)>=6; print $_"""
  conda:
    "env/align.yml"
  shell:
    "ml samtools; samtools view -b -F 0X800 -F 0X400 {input} |samtools collate -O - | bedtools bamtobed -bedpe -i stdin| perl -lane {params.cmd:q} > {output}"

rule markdup:
  input:
    "bam/{sra}_fastp.sort.bam"
  output:
    bam = "bam/{sra}_dedup.bam",
    m = "bam/{sra}_picard.txt"
  conda:
    "env/align.yml"
  shell:
    "ml picard; java -jar {PICARD} MarkDuplicates I={input} O={output.bam} M={output.m}"

# rule sort_bam:
#   input:
#     "bam/{sra}_fastp.bam"
#   output:
#     "bam/{sra}_fastp.sort.bam"
#   shell:
#     "ml samtools;samtools sort -o {output} {input}"

rule bwa:
  input:
    expand("fastp_fq/{sra}_fastp_{n}.fastq", n=[1, 2], allow_missing=True)
  output:
    "bam/{sra}_fastp.sort.bam"
  threads:
    4
  conda:
    "env/align.yml"
  shell:
    """ml bwa samtools;bwa mem -5SP -t {threads} {HG38_MM10} {input[0]} {input[1]} | samtools sort -o {output}"""

rule fastp:
  input:
    expand("fastq/{sra}_mosaic_dust_{n}.fastq", n=[1, 2], allow_missing=True)
  output:
    fq = expand("fastp_fq/{sra}_fastp_{n}.fastq", n=[1, 2], allow_missing=True),
    json = "fastp_fq/{sra}_fastp.json",
    html = "fastp_fq/{sra}_fastp.html"
  conda:
    "env/align.yml"
  shell:
    """fastp --dedup -i {input[0]} -I {input[1]} -j {output.json} -h {output.html} -o {output.fq[0]} -O {output.fq[1]}"""
