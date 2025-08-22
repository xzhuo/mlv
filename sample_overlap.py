# Create a sample overlap matrix to help decide which samples should be combined in analysis.

import itertools
import os
import argparse
import re

"""
Sample input lines:
chr1	13746	13892	SRR13105370.14195290
chr1	200525	200661	SRR13105360.17618327,SRR14032981.17618327
chr1	273015	273112	SRR13105370.14268687
chr1	774282	774372	ERR1856018.63313094
chr1	778338	778525	SRR16771772.4863824,SRR29662683.33948916
chr1	810935	811070	SRR13105375.51822124
"""

def process_bed(input, out):
    matrix = {} # the 2d matrix
    all_samples = set()
    with open(input, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split()
            if len(fields) < 4:
                continue
            reads = fields[3]
            read_list = reads.split(',')
            sample_list = map(lambda x: x.split('.')[0], read_list)
            sample_set = set(sample_list)
            all_samples.update(sample_set)
            if len(sample_set) > 1:
                for sample_pair in itertools.combinations(sample_set, 2):
                    # breakpoint()
                    paired_list = list(sample_pair)
                    paired_list.sort()
                    try:
                        matrix[paired_list[0]][paired_list[1]] += 1
                    except KeyError:
                        matrix[paired_list[0]] = {paired_list[1]: 1}
            else:
                one = sample_set.pop()
                try:
                    matrix[one]["NA"] +=1
                except KeyError:
                    matrix[one] = {"NA": 1}
    for sample in all_samples:
        if sample not in matrix:
            matrix[sample] = {}
        for other_sample in all_samples:
            if other_sample not in matrix[sample]:
                matrix[sample][other_sample] = 0
            if "NA" not in matrix[sample]:
                matrix[sample]["NA"] = 0

    with open(out, 'w') as f:
        f.write("Sample\t" + "\t".join(str(x) for x in matrix.keys()) + "\tNA"+ "\n")
        for row in matrix:
            line_list = list(map(lambda x: matrix[row][x], matrix.keys()))
            line_list.append(matrix[row]["NA"])
            f.write(str(row) + "\t" + "\t".join(str(x) for x in line_list) + "\n")


def main():
    parser = argparse.ArgumentParser(description='parse a bedpe file and summerize insertion strands')
    parser.add_argument('-i', '--input', type=str, required=True,
                        help='input merged bed file, with list ofreads as comma separated values in column 4')
    parser.add_argument('-o', '--out', type=str, required=True,
                        help='output summery txt file')

    args = parser.parse_args()
    bedpe_file = os.path.abspath(args.input)
    if not os.path.exists(bedpe_file):
        raise ValueError("--input file does not exist!")
    process_bed(args.input, args.out)

if __name__ == '__main__':
    main()