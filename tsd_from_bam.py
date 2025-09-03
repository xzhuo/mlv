import os
import argparse
import pysam
from Bio.Seq import Seq
from parse_insertion_bedpe import Read


def extract_tsd(bam_file, input_file, out):
    bam_file = pysam.AlignmentFile(bam_file, "rb")
    all_reads = pysam.IndexedReads(bam_file, False) # close it since I am not gonna change the bam.
    all_reads.build()
    for read in all_reads.find("SRR13105356.11541486"):
        if read.is_unmapped or read.is_secondary or read.is_duplicate or read.is_supplementary:
            continue
        breakpoint()

    with open(input_file, 'r') as f:
        for line in f:
            if line.startswith('#') or line.strip() == "":
                continue
            bedpe_line = Read(line)
            fetch_human = bedpe_line.soft_clip_1
            fetch_mouse = bedpe_line.soft_clip_2
            human_cigar = bedpe_line.cigar_1
            mouse_cigar = bedpe_line.cigar_2
            read_id = bedpe_line.name
            if bedpe_line.position == "up":
                if bedpe_line.soft_clip_1:
                    cigar = bedpe_line.cigar_1
                    boundary_seq = Seq(read.query_sequence)[read.query_alignment_end - 10:read.query_alignment_end + 10]
                if bedpe_line.soft_clip_2:
                    cigar = bedpe_line.cigar_2
                    full_read = Seq(read.get_forward_sequence()).reverse_complement()
            elif bedpe_line.position == "down":
                if bedpe_line.soft_clip_1:
                    cigar = bedpe_line.cigar_1
                    full_read = Seq(read.query_sequence)[read.query_alignment_start - 10:read.query_alignment_start + 10]
                if bedpe_line.soft_clip_2:
                    cigar = bedpe_line.cigar_2
                    full_read = Seq(read.get_forward_sequence())
            else:
                raise ValueError("wrong strand info in the bedpe file!")
            # # up or down
            # # Check for soft clipping
            # if read.cigar[0][0] == 4:  # Soft clip
        #     tsd = read.query_sequence[read.cigar[0][1]:]
        #     fout.write(f">{read.query_name}\n{tsd}\n")

def main():
    parser = argparse.ArgumentParser(description='Extract tsd sequences from a bam file')
    parser.add_argument('-b', '--bam', type=str, required=True,
                        help='input bam file')
    parser.add_argument('-i', '--input', type=str, required=False,
                        help='input bedpe file')
    parser.add_argument('-o', '--out', type=str, required=False,
                        help='output extracted seq')


    args = parser.parse_args()
    input_file = os.path.abspath(args.input)
    bam_file = os.path.abspath(args.bam)
    # if not os.path.exists(bam_file) or not os.path.exists(input_file):
    #     raise ValueError("--bam file or --input file does not exist!")
    extract_tsd(bam_file, input_file, args.out)

if __name__ == '__main__':
    main()
