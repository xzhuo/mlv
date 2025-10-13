import os
import argparse
import pysam
import re

def parse_fa(all_bam_reads, input_file, output_file, failed_file):
    with open(input_file, 'r') as in_put:
        with open(output_file, 'w') as out:
            with open(failed_file, 'w') as failed:
                insertion = [] # chrom, start, end, strand, read_num, inconsitend_read_num, TSD_annotation
                read = [] # read ID, strand, genome
                for line in in_put:
                    line = line.strip()
                    if line == "":
                        continue
                    elif line.startswith('#'):
                        if insertion:
                            insertion[1] = str(sorted(insertion_pos["+"].items(), key=lambda item: item[1], reverse=True))
                            insertion[2] = str(sorted(insertion_pos["-"].items(), key=lambda item: item[1], reverse=True))
                            out.write('\t'.join(insertion) + '\n')
                        line = line.removeprefix('# ')
                        insertion = line.split()
                        insertion_pos = {"+":{}, "-":{}}
                        # out.write(line + '\n')
                        failed.write(line + '\n')
                    elif line.startswith('>'):
                        line = line.removeprefix('>')
                        read = line.split('_') # SRR21439924.8076885_+_mm10
                    else:
                        matches = list(re.finditer(r'\|', line))
                        # move on and print the line to a new file with no matches;
                        if len(matches) == 0:
                            failed.write(line + '\n')
                        # raise error if there are more than 1 matches;
                        elif len(matches) > 1:
                            raise ValueError("more than 1 matches found in the line: " + line)
                        elif len(matches) == 1:
                            # raise error if "TTTCA|" in a + strand match;
                            if read[2] == '+' and "TTTCA|" in line:
                                raise ValueError("TTTCA| found in a + strand match: " + line)
                            # raise error if "|TGAAA" in a - strand match;
                            elif read[2] == '-' and "|TGAAA" in line:
                                raise ValueError("|TGAAA found in a - strand match: " + line)
                            else:
                                match = matches.pop()
                                offset = match.start() - 20
                                read_id, strand, genome = read
                                tsd_pos = None
                                if genome =="hg38":
                                    bam_reads = [read for read in all_bam_reads.find(read_id) if re.match("GRCh38",read.reference_name) and not (read.is_unmapped or read.is_secondary or read.is_duplicate or read.is_supplementary)]
                                    human_read = bam_reads.pop() # I don't need to check if there are multiple reads since it has been checked in the previous script.
                                    if strand == "+":
                                        ori_pos = human_read.reference_end
                                        tsd_pos = ori_pos + offset
                                    elif strand == "-":
                                        ori_pos = human_read.reference_start
                                        tsd_pos = ori_pos - offset

                                elif genome == "mm10":
                                    bam_reads = [read for read in all_bam_reads.find(read_id) if re.match("mm10",read.reference_name) and not (read.is_unmapped or read.is_secondary or read.is_duplicate or read.is_supplementary)]
                                    mouse_read = bam_reads.pop() # I don't need to check if there are multiple reads since it has been checked in the previous script
                                    supp_reads = [read for read in all_bam_reads.find(read_id) if re.match("GRCh38",read.reference_name) and read.is_supplementary and read.is_read1 == mouse_read.is_read1]
                                    if supp_reads:
                                        supp_read = supp_reads.pop()
                                        human_clip_length = mouse_read.query_alignment_start if mouse_read.is_reverse else mouse_read.query_length - mouse_read.query_alignment_end
                                        if strand == "+":
                                            human_clip_length += offset
                                            supp_human_pos = supp_read.reference_end
                                            tsd_pos = supp_human_pos - (human_clip_length - (supp_read.query_length - supp_read.query_alignment_start))
                                        elif strand == "-":
                                            human_clip_length -= offset
                                            supp_human_pos = supp_read.reference_start
                                            tsd_pos = supp_human_pos + (human_clip_length - supp_read.query_alignment_end)
                                    else:
                                        tsd_pos = "?"
                                if tsd_pos:
                                    try:
                                        insertion_pos[read[1]][str(tsd_pos)] += 1
                                    except KeyError:
                                        insertion_pos[read[1]][str(tsd_pos)] = 1
                                    # out.write('\t'.join([read_id, strand, genome, str(tsd_pos)]) + '\n')
                        else:
                            raise ValueError("unexpected error in the line: " + line)
    
                if insertion:
                    insertion[1] = str(sorted(insertion_pos["+"].items(), key=lambda item: item[1], reverse=True))
                    insertion[2] = str(sorted(insertion_pos["-"].items(), key=lambda item: item[1], reverse=True))
                    out.write('\t'.join(insertion) + '\n')

def main():
    parser = argparse.ArgumentParser(description='summarize the manual check results in a fasta file to a table')
    parser.add_argument('-b', '--bam', type=str, required=True,
                        help='input bam file')
    parser.add_argument('-i', '--input', type=str, required=False,
                        help='input manually checked fasta file')
    parser.add_argument('-o', '--out', type=str, required=False,
                        help='output summary of TSD annotations')
    parser.add_argument('-f', '--failed', type=str, required=False,
                        help='output fasta of reads with out TSD boundary')

    args = parser.parse_args()
    input_file = os.path.abspath(args.input)
    bam_file = os.path.abspath(args.bam)
    # if not os.path.exists(bam_file) or not os.path.exists(input_file):
    #     raise ValueError("--bam file or --input file does not exist!")
    bam_file = pysam.AlignmentFile(bam_file, "rb")
    all_bam_reads = pysam.IndexedReads(bam_file, False) # close it since I am not gonna change the bam.
    all_bam_reads.build()
    parse_fa(all_bam_reads, input_file, args.out, args.failed)

if __name__ == '__main__':
    main()
