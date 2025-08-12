import os
import argparse

def summarize_insertions(bedpe_file, out, failed, distance):
    '''
    five lines of bedpe file, must be sorted:
    chr1	10105675	10105776	RLTR_Mm	0	84	DRR349015.15272214	+	-
    chr1	10105675	10105776	RLTR_Mm	0	84	DRR349015.9109586	+	-
    chr1	10105701	10105785	RLTR_Mm	0	59	DRR349015.17456506	+	-
    chr14	20454045	20454136	RLTR_Mm	0	54	DRR349014.4907801	+	-
    chr14	20454124	20454225	RLTR_Mm	686	742	DRR349012.20595080	-	+
    '''

    curr_chr = ''
    curr_pos = int()
    insertion_strand = ''
    curr_strand = ''
    insertion_number = int()
    consistent_read = int()
    inconsistent_read = int()
    failed_list = []
    curr_num_reads = int()
    with open(bedpe_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split()
            if len(fields) < 9:
                continue
            chrom_1, start_1, end_1, chrom_2, start_2, end_2, name, strand_1, strand_2 = fields[:9]
            start_1, end_1,start_2, end_2= map(int, [start_1, end_1, start_2, end_2])
            if chrom_1 == curr_chr and start_1 < int(curr_pos) + distance:
                # existing insertion
                curr_pos = end_1
                line_insertion = "-" if strand_1 == strand_2 else '+'
                curr_num_reads += 1
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
                curr_chr = chrom_1
                curr_pos = end_1
                insertion_strand = "-" if strand_1 == strand_2 else '+'
                curr_strand = strand_1
                insertion_number += 1 if curr_num_reads > 1 else 0
                curr_num_reads = 1
        insertion_number += 1 if curr_num_reads > 1 else 0

    with open(out, 'w') as f:
        f.write(f"{insertion_number}\t{consistent_read}\t{inconsistent_read}\n")
    if failed:
        with open(failed, 'w') as failed_out:
            for failed_line in failed_list:
                failed_out.write(f"{failed_line}\n")


def main():
    parser = argparse.ArgumentParser(description='parse a bedpe file and summerize insertion strands')
    parser.add_argument('-i', '--input', type=str, required=True,
                        help='input bedpe file')
    parser.add_argument('-o', '--out', type=str, required=True,
                        help='output summery txt file')
    parser.add_argument('-f', '--failed', type=str, required=False,
                        help='output seq does not match')
    parser.add_argument('-d', '--distance', type=int, default=1000,
                        help='distance within which we consider the same insertion')

    args = parser.parse_args()
    bedpe_file = os.path.abspath(args.input)
    if not os.path.exists(bedpe_file):
        raise ValueError("--input file does not exist!")
    summarize_insertions(bedpe_file, args.out, args.failed, args.distance)

if __name__ == '__main__':
    main()