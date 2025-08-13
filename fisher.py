import os
import argparse
import pybedtools
from itertools import repeat
from multiprocessing import Pool

def fisher_test(input, te_bed, genomesize_file, sample, te):
    print("working on: {} in sample {}".format(te, sample))

    # filter the input bed: only keep the rows that match the specified sample
    sample_bed = input.filter(lambda x: x[3].startswith(sample))
    # fisher exact test
    fisher_output = sample_bed.fisher(te_bed, g=genomesize_file)
    ratio = fisher_output.ratio
    two_tail = fisher_output.two_tail
    right_tail = fisher_output.right_tail
    left_tail = fisher_output.left_tail
    intersect_number = fisher_output.table['in -a']['in -b']
    return [sample, te, intersect_number, ratio, two_tail, left_tail, right_tail]
    # out_line = "\t".join(str(x) for x in [sample, te, intersect_number, ratio, two_tail, left_tail, right_tail])
    # return out_line


def main():
    parser = argparse.ArgumentParser(description='run multiple bedtools fisher exact test')
    parser.add_argument('-m', '--mei', type=str, default="RLTR4_Mm",
                        help='The mobile element to test')
    parser.add_argument('-i', '--input', type=str, required=True,
                        help='input bed as a for bedtools.')
    parser.add_argument('-r', '--rmsk', type=str, required=True,
                        help='rmsk bed input as b for bedtools')
    parser.add_argument('-g', '--genomesize', type=str, required=True,
                        help='genome size file')
    parser.add_argument('-o', '--out', type=str, required=True,
                        help='output fisher test results')
    parser.add_argument('-l', '--list', type=str, required=True,
                        help='output a list of all reads that intersect with the TE')
    parser.add_argument('-t', '--threads', type=int, default=1,
                        help='multi-threading')

    args = parser.parse_args()

    input_file = os.path.abspath(args.input)
    rmsk_file = os.path.abspath(args.rmsk)
    genomesize_file = os.path.abspath(args.genomesize)

    if not os.path.exists(input_file):
        raise ValueError("--input file does not exist!")
    if not os.path.exists(rmsk_file):
        raise ValueError("--rmsk file does not exist!")
    if not os.path.exists(genomesize_file):
        raise ValueError("--genome size file does not exist!")

    te = args.mei

    # print the list of samples that intersect with the specified mobile element.
    input = pybedtools.BedTool(input_file)
    rmsk = pybedtools.BedTool(rmsk_file)
    # filter rmsk: only keep the rows that match the specified mobile element
    te_bed = rmsk.filter(lambda x: x[3] == te)
    # run bedtools intersect to get the intersection of the input bed and the TE bed
    intersect_bed = input.intersect(te_bed, u=True)
    read_ids = [str(x.name) for x in intersect_bed]
    with open(args.list, 'w') as lst:
        lst.write("\n".join(read_ids) + "\n")

    # parse the input file and generate a unique list of samples
    sample_list = set()
    with open(input_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 4:
                continue
            sample, id = fields[3].split('.', 2)  # assuming the sample with id is in the 4th column separated by a dot
            sample_list.add(sample)
    print("{} samples found in input file".format(len(sample_list)))
    results = []

    if args.threads == 1:
        for sample in sample_list:
            out_line = fisher_test(input, te_bed, genomesize_file, sample)
            results.append(out_line)
    else:
        with Pool(processes=args.threads) as pool:
            results = pool.starmap(fisher_test, map(lambda x:(input, te_bed, genomesize_file, x),sample_list))

    with open(args.out, 'w') as out:
        for each_out in results:
            out_line = "\t".join(str(x) for x in each_out)  # format the output line
            out.write(out_line + "\n")

if __name__ == '__main__':
    main()
