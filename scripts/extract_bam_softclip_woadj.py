import os
import argparse
import pysam
import re
# from Bio.Seq import Seq


def extract_softclip(bam_file, out, part):
    bam = pysam.AlignmentFile(bam_file, "r")
    with open(out, 'w') as outfile:
        # Iterate through each read in the BAM file 
        for read in bam:
            if read.is_supplementary or read.is_duplicate or read.is_secondary or read.is_unmapped:
                continue
            else:
                if (read.is_reverse and read.query_alignment_start > 0) or (read.is_forward and read.query_alignment_end < read.query_length):
                    offset = read.query_length - read.query_alignment_start if read.is_reverse else read.query_alignment_end
                    forward_sequence = read.get_forward_sequence()
                    subseq = forward_sequence[offset:] if part == 'clip' else forward_sequence[:offset]
                    sa = read.get_tag('SA') if read.has_tag('SA') else None
                    outfile.write(f">{read.query_name}_clip\t# {sa}\n{subseq}\n")


def main():
    parser = argparse.ArgumentParser(description='Extract soft clipped sequences from all reads of a bam file')
    parser.add_argument('-b', '--bam', type=str, required=True,
                        help='input bam file')
    parser.add_argument('-o', '--out', type=str, required=True,
                        help='output extracted all softclipped seq in fasta format')
    parser.add_argument('-p', '--part', type=str, choices=["clip", "align"], required=True,
                        help='extract softclipped seq or the subsequence with softclipped region removed, options: clip or align')

    args = parser.parse_args()
    bam_file = os.path.abspath(args.bam)
    if not os.path.exists(bam_file):
        raise ValueError("--bam file does not exist!")
    if args.part not in ['clip', 'align']:
        raise ValueError("--part must be either 'clip' or 'align'")
    extract_softclip(bam_file, args.out, args.part)

if __name__ == '__main__':
    main()
