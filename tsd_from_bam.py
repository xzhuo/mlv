import os
import argparse
import pysam
import re
from Bio.Seq import Seq
from parse_insertion_bedpe import Read

def boundary_subseq(read, direction, revcom = False, distance = 10): # direct has to be up or down.
    boundary_subseq = ""
    if direction == "up":
        boundary_start = read.query_alignment_end - distance if read.query_alignment_end - distance > 0 else 0
        gap_start = "" if read.query_alignment_end - distance > 0 else "-" * (distance - read.query_alignment_end)
        boundary_end = read.query_alignment_end + distance if read.query_alignment_end + distance < read.query_length else read.query_length
        gap_end = "" if read.query_alignment_end + distance < read.query_length else "-" * (distance - (read.query_length - read.query_alignment_end))
        boundary_subseq = Seq(gap_start + read.query_sequence[boundary_start:boundary_end] + gap_end).reverse_complement() if revcom else Seq(gap_start + read.query_sequence[boundary_start:boundary_end] + gap_end)
    elif direction == "down":
        boundary_start = read.query_alignment_start - distance if read.query_alignment_start - distance > 0 else 0
        gap_start = "" if read.query_alignment_start - distance > 0 else "-" * (distance - read.query_alignment_start)
        boundary_end = read.query_alignment_start + distance if read.query_alignment_start + distance < read.query_length else read.query_length
        gap_end = "" if read.query_alignment_start + distance < read.query_length else "-" * (distance - (read.query_length - read.query_alignment_start))
        boundary_subseq = Seq(gap_start + read.query_sequence[boundary_start:boundary_end] + gap_end).reverse_complement() if revcom else Seq(gap_start + read.query_sequence[boundary_start:boundary_end] + gap_end)
    else:
        raise ValueError("Invalid direction")

    return boundary_subseq

def extract_tsd(bam_file, input_file, out):
    bam_file = pysam.AlignmentFile(bam_file, "rb")
    all_reads = pysam.IndexedReads(bam_file, False) # close it since I am not gonna change the bam.
    all_reads.build()

    with open(input_file, 'r') as f:
        with open(out, 'w') as out:
            for line in f:
                if line.startswith('#'):
                    out.write(line)
                    continue
                if line.strip() == "":
                    continue

                bedpe_line = Read(line)
                read_id = bedpe_line.name
                human_read_list = [read for read in all_reads.find(read_id) if re.match("GRCh38",read.reference_name) and not (read.is_unmapped or read.is_secondary or read.is_duplicate or read.is_supplementary)]
                mouse_read_list = [read for read in all_reads.find(read_id) if re.match("mm10",read.reference_name) and not (read.is_unmapped or read.is_secondary or read.is_duplicate or read.is_supplementary)]
                # There should be only one read in each list.
                try:
                    human_read = human_read_list.pop() # an pysam.AlignedSegment object
                    mouse_read = mouse_read_list.pop() # an pysam.AlignedSegment object
                except IndexError:
                    raise IndexError(f"cannot find the proper paired read {read_id} in the bam file!")
                if human_read_list or mouse_read_list:
                    raise ValueError(f"find more than a paired reads {read_id} in the bam file!")
                fast_list = []
                if bedpe_line.position == "up":
                    if bedpe_line.soft_clip_1 and re.search(r'(\d+)S$', bedpe_line.cigar_1):
                        boundary_seq = boundary_subseq(human_read, "up", False, 10)
                        fasta = f">{read_id}_+_hg38\n{boundary_seq}"
                        fast_list.append(fasta)
                    if bedpe_line.soft_clip_2:
                        if (mouse_read.is_reverse and re.search(r'^(\d+)S', bedpe_line.cigar_2)) or (mouse_read.is_forward and re.search(r'(\d+)S$', bedpe_line.cigar_2)):
                            direction = "down" if mouse_read.is_reverse else "up"
                            revcom = False if mouse_read.is_reverse else True
                            boundary_seq = boundary_subseq(mouse_read, direction, revcom, 10)
                            fasta = f">{read_id}_+_mm10\n{boundary_seq}"
                            fast_list.append(fasta)
                elif bedpe_line.position == "down":
                    if bedpe_line.soft_clip_1 and re.search(r'^(\d+)S', bedpe_line.cigar_1):
                        boundary_seq = boundary_subseq(human_read, "down", False, 10)
                        fasta = f">{read_id}_-_hg38\n{boundary_seq}"
                        fast_list.append(fasta)
                    if bedpe_line.soft_clip_2:
                        if (mouse_read.is_reverse and re.search(r'^(\d+)S', bedpe_line.cigar_2)) or (mouse_read.is_forward and re.search(r'(\d+)S$', bedpe_line.cigar_2)):
                            direction = "down" if mouse_read.is_reverse else "up"
                            revcom = True if mouse_read.is_reverse else False
                            boundary_seq = boundary_subseq(mouse_read, direction, revcom, 10)
                            fasta = f">{read_id}_-_mm10\n{boundary_seq}"
                            fast_list.append(fasta)
                else:
                    raise ValueError("wrong strand info in the bedpe file!")
                if fast_list:
                    out.write("\n".join(fast_list) + "\n")

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
