import os
import argparse
import re

class Read:
    def __init__(self, line):
        self.line = line
        self.dissect_line()

    def dissect_line(self):
        fields = self.line.strip().split()
        self.chrom_1, start_1, end_1, self.chrom_2, start_2, end_2, self.name, self.strand_1, self.strand_2, self.cigar_1, self.cigar_2 = fields[:11]
        self.start_1, self.end_1, self.start_2, self.end_2= map(int, [start_1, end_1, start_2, end_2])
        self.insertion_strand = "-" if self.strand_1 == self.strand_2 else '+'
        if self.strand_1 == '+':
            self.position = "up"
            self.upstream_boundary = self.end_1 
        if self.strand_1 == '-':
            self.position = "down"
            self.downstream_boundary = self.start_1
        if re.search(r'(\d+)S', self.cigar_1):
            self.soft_clip_1 = True
        if re.search(r'(\d+)S', self.cigar_2):
            self.soft_clip_2 = True
        self.any_soft_clip = self.soft_clip_1 or self.soft_clip_2




class Insertion:
    def __init__(self):
        self.chrom = ''
        self.strand = ''
        self.insertion_strand = ''
        self.up_reads = []
        self.down_reads = []
        self.upstream_boundary = int()
        self.downstream_boundary = int()
        self.soft_clipped_reads = []
    
    def add_read(self, read):
        self.chrom = read.chrom_1
        self.strand = read.strand_1
        self.insertion_strand = read.insertion_strand
        
        if read.position == "up":
            self.up_reads.append(read)
            self.up_softclip = self.up_softclip and read.soft_clip_1
            if read.upstream_boundary > self.upstream_boundary:
                self.upstream_boundary = read.upstream_boundary
        
        if read.position == "down":
            self.down_reads.append(read)
            self.down_softclip = self.down_softclip and read.soft_clip_1
            if read.downstream_boundary < self.downstream_boundary:
                self.downstream_boundary = read.downstream_boundary
    
    def does_contain(self, read, distance):
        if self.chrom == read.chrom_1 and read.start_1 < self.last_pos + distance:
            pass
            if self.insertion_strand == read.insertion_strand:
                pass
        else:
            return False


def summarize_insertions(bedpe_file, out, failed, tsd, distance):
    '''
    five lines of bedpe file, must be sorted:
    chr1	1014478	1014568	emv2_LTR	183	273	ERR1856015.51879930	+	+	90M	90M
    chr1	1779603	1779693	emv2_LTR	83	173	ERR1856016.58729964	-	-	90M	90M
    chr1	2411502	2411592	emv2_LTR	422	512	ERR1856015.150948519	+	+	90M	90M
    chr1	9673014	9673104	emv2_LTR	383	473	ERR1856017.131318684	+	+	90M	90M
    '''

    insertions = [] # a list of Insertion class.
    with open(bedpe_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue

            read = Read(line)
            if len(insertions) > 0 and insertions[-1].does_contain(read, distance):
                insertions[-1].add_read(read)
            else:
                insertion = Insertion()
                insertion.add_read(read)
                insertions.append(insertion)


            if chrom_1 == curr_chr and start_1 < int(curr_pos) + distance:
                # existing insertion
                curr_pos = end_1
                line_insertion = "-" if strand_1 == strand_2 else '+'
                curr_num_reads += 1
                if strand_1 == '+' and (re.search(r'(\d+)S', cigar_1) or re.search(r'(\d+)S', cigar_2)):
                    curr_left_softclip.append(line.strip())
                if strand_1 == '-' and (re.search(r'(\d+)S', cigar_1) or re.search(r'(\d+)S', cigar_2)):
                    curr_right_softclip.append(line.strip())
                if (curr_strand == "-" and strand_1 == '+') or (line_insertion != insertion_strand):
                    inconsistent_read += 1
                    failed_list.append(line.strip())
                    # next try to report inconsistent reads per insertion.
                else:
                    insertion_strand = line_insertion
                    curr_strand = strand_1
                    consistent_read += 1
            else:
                # new insertion
                if len(curr_left_softclip) > 0 and len(curr_right_softclip) > 0:
                    # report soft-clipped reads
                    tsd_list.extend(curr_left_softclip + curr_right_softclip)
                curr_left_softclip = []
                curr_right_softclip = []
                if strand_1 == '+' and (re.search(r'(\d+)S', cigar_1) or re.search(r'(\d+)S', cigar_2)):
                    curr_left_softclip.append(line.strip())
                if strand_1 == '-' and (re.search(r'(\d+)S', cigar_1) or re.search(r'(\d+)S', cigar_2)):
                    curr_right_softclip.append(line.strip())
                curr_chr = chrom_1
                curr_pos = end_1
                insertion_strand = "-" if strand_1 == strand_2 else '+'
                curr_strand = strand_1
                insertion_number += 1 if curr_num_reads > 1 else 0
                curr_num_reads = 1

        if len(curr_left_softclip) > 0 and len(curr_right_softclip) > 0:
            # report soft-clipped reads
            tsd_list.append(curr_left_softclip + curr_right_softclip)
        insertion_number += 1 if curr_num_reads > 1 else 0

    with open(out, 'w') as f:
        f.write(f"{insertion_number}\t{consistent_read}\t{inconsistent_read}\n")
    if failed:
        with open(failed, 'w') as failed_out:
            for failed_line in failed_list:
                failed_out.write(f"{failed_line}\n")
    if tsd:
        with open(tsd, 'w') as tsd_out:
            for tsd_line in tsd_list:
                tsd_out.write(f"{tsd_line}\n")


def main():
    parser = argparse.ArgumentParser(description='parse a bedpe file and summerize insertion strands')
    parser.add_argument('-i', '--input', type=str, required=True,
                        help='input bedpe file')
    parser.add_argument('-o', '--out', type=str, required=True,
                        help='output summery txt file')
    parser.add_argument('-f', '--failed', type=str, required=False,
                        help='output seq does not match')
    parser.add_argument('-t', '--tsd', type=str, required=False,
                        help='output TSD seq')
    parser.add_argument('-d', '--distance', type=int, default=1000,
                        help='distance within which we consider the same insertion')

    args = parser.parse_args()
    bedpe_file = os.path.abspath(args.input)
    if not os.path.exists(bedpe_file):
        raise ValueError("--input file does not exist!")
    summarize_insertions(bedpe_file, args.out, args.failed, args.tsd, args.distance)

if __name__ == '__main__':
    main()