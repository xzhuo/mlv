import os
import argparse
import re
# from statistics import median
# from functools import reduce

class Read:
    def __init__(self, line):
        self.line = line.strip()
        self.dissect_line()
        self.edge_score = int()
        self.closest = False
        self.at_edge = False
        self.inconsistent = False

    def dissect_line(self):
        fields = self.line.split()
        self.chrom_1, start_1, end_1, self.chrom_2, start_2, end_2, self.name, self.strand_1, self.strand_2, self.cigar_1, self.cigar_2 = fields[:11]
        self.start_1, self.end_1, self.start_2, self.end_2= map(int, [start_1, end_1, start_2, end_2])
        self.insertion_strand = "-" if self.strand_1 == self.strand_2 else '+'
        if self.strand_1 == '+':
            self.position = "up"
            self.up_edge = self.end_1 
        if self.strand_1 == '-':
            self.position = "down"
            self.down_edge = self.start_1
        self.soft_clip_1 = True if re.search(r'(\d+)S', self.cigar_1) else False
        self.soft_clip_2 = True if re.search(r'(\d+)S', self.cigar_2) else False
        self.any_soft_clip = self.soft_clip_1 or self.soft_clip_2

    def print(self):
        return self.line

class Insertion:
    def __init__(self):
        self.chrom = ''
        self.insertion_strand = ''
        self.up_reads = []
        self.down_reads = []
        self.last_pos = int()
        self.at_up_edge = False
        self.at_down_edge = False
        self.up_tsd = False
        self.down_tsd = False
        self.max_up_score = int()
        self.max_down_score = int()
        self.at_up_edge = False
        self.at_down_edge = False
        self.up_closest = int()
        self.down_closest = int()


    def add_read(self, read):
        self.chrom = read.chrom_1
        self.insertion_strand = read.insertion_strand
        self.last_pos = max(self.last_pos, read.end_1)
        
        if read.position == "up":
            self.up_reads.append(read)
            # if read.soft_clip_1:
            #     self.up_edgies_softclip.append(read.up_edge)
            # self.up_softclip = self.up_softclip and read.soft_clip_1
            # self.up_limit = max(self.up_limit, read.up_edge)
            # if read.soft_clip_1:
            #     self.up_edgies.append(read.up_edge)
            if read.any_soft_clip:
                self.up_tsd = True

        if read.position == "down":
            self.down_reads.append(read)
            # self.down_softclip = self.down_softclip and read.soft_clip_1
            # if self.down_limit == 0:
            #     self.down_limit = read.down_edge
            # else:
            #     self.down_limit = min(self.down_limit, read.down_edge)
            # if read.soft_clip_1:
            #     self.down_edgies.append(read.down_edge)
            if read.any_soft_clip:
                self.down_tsd = True

    def resolve_insertion_site(self, tsd = 4, wiggle_distance = 20):
        for read in self.up_reads:
            conflict_up_num = len([i for i in self.up_reads if i.up_edge > read.up_edge + wiggle_distance])
            conflict_down_num = len([i for i in self.down_reads if i.down_edge < read.up_edge - tsd - wiggle_distance])
            conflict_num = conflict_up_num + conflict_down_num
            read.edge_score = 100 / (conflict_num + 1) # Convert the number of other reads conflict with it being at the edge to a score. 100 / (conflict_num + 1)

        for read in self.down_reads:
            conflict_up_num = len([i for i in self.up_reads if i.up_edge > read.down_edge + tsd + wiggle_distance])
            conflict_down_num = len([i for i in self.down_reads if i.down_edge < read.down_edge - wiggle_distance])
            conflict_num = conflict_up_num + conflict_down_num
            read.edge_score = 100 / (conflict_num + 1) # Convert the number of other reads conflict with it being at the edge to a score. 100 / (conflict_num + 1)
        if len(self.up_reads) > 0:
            self.max_up_score = max([read.edge_score for read in self.up_reads])
        if len(self.down_reads) > 0:
            self.max_down_score = max([read.edge_score for read in self.down_reads])

        # label closest, at_edge, and inconsistent on all reads,
        # label upstream and downstream closest position and whether the edge is soft clipped on the insertion.
        up_closest_all = []
        down_closest_all = []
        for read in self.up_reads:
            if read.edge_score == self.max_up_score:
                read.closest = True
                up_closest_all.append(read.up_edge)
                if read.soft_clip_1:
                    read.at_edge = True
                    self.at_up_edge = True
        for read in self.down_reads:
            if read.edge_score == self.max_down_score:
                read.closest = True
                down_closest_all.append(read.down_edge)
                for i in self.up_reads:
                    if i.up_edge > read.down_edge + tsd + wiggle_distance:
                        i.inconsistent = True
                for i in self.down_reads:
                    if i.down_edge < read.down_edge - wiggle_distance:
                        i.inconsistent = True
                if read.soft_clip_1:
                    read.at_edge = True
                    self.at_down_edge = True
        self.up_closest = max(up_closest_all) if len(up_closest_all) > 0 else 0
        self.down_closest = min(down_closest_all) if len(down_closest_all) > 0 else 0
        for read in self.up_reads:
            if read.up_edge > self.up_closest + wiggle_distance:
                read.inconsistent = True
        for read in self.down_reads:
            if read.down_edge < self.down_closest - wiggle_distance:
                read.inconsistent = True

    def reset_scores(self):
        self.max_up_score = 0
        self.max_down_score = 0
        self.min_up_score = 0
        self.min_down_score = 0
        self.up_closest = 0
        self.down_closest = 0
        for read in self.up_reads:
            read.edge_score = 0
            read.closest = False
            read.at_edge = False
            read.inconsistent = False
        for read in self.down_reads:
            read.edge_score = 0
            read.closest = False
            read.at_edge = False
            read.inconsistent = False

    def diagnose(self, wiggle_distance = 20):
        up_softclip_pos = [read.up_edge for read in self.up_reads if read.soft_clip_1]
        down_softclip_pos = [read.down_edge for read in self.down_reads if read.soft_clip_1]
        return up_softclip_pos, down_softclip_pos

    def find_secondary_up_insertion(self, tsd = 4, wiggle_distance = 20):
        if len(self.get_inconsistent_up_reads()) > 2:
            # If there are more than 2 inconsistent reads, we can try to find a secondary insertion
            new_insertion = Insertion()
            for read in self.get_inconsistent_up_reads():
                try:
                    self.up_reads.remove(read)
                except ValueError:
                    raise ValueError("Downstream read not found in the list up_reads")
                new_insertion.add_read(read)
            # find all reads from the old insertion that is compatible with the new one
            new_insertion.resolve_insertion_site()
            separate_reads = [read for read in self.down_reads if read.down_edge > new_insertion.up_closest - tsd - wiggle_distance]
            for read in separate_reads:
                try:
                    self.down_reads.remove(read)
                except ValueError:
                    raise ValueError("Downstream read not found in the list down_reads")
                new_insertion.add_read(read)
            new_insertion.resolve_insertion_site()
            self.resolve_insertion_site()
            return new_insertion

    def find_secondary_down_insertion(self, tsd = 4, wiggle_distance = 20):
        if len(self.get_inconsistent_down_reads()) > 2:
            # If there are more than 2 inconsistent reads, we can try to find a secondary insertion
            new_insertion = Insertion()
            for read in self.get_inconsistent_down_reads():
                try:
                    self.down_reads.remove(read)
                except ValueError:
                    raise ValueError("Downstream read not found in the list down_reads")
                new_insertion.add_read(read)
            # find all reads from the old insertion that is compatible with the new one
            new_insertion.resolve_insertion_site()
            separate_reads = [read for read in self.up_reads if read.up_edge < new_insertion.down_closest + tsd + wiggle_distance]
            for read in separate_reads:
                try:
                    self.up_reads.remove(read)
                except ValueError:
                    raise ValueError("Upstream read not found in the list up_reads")
                new_insertion.add_read(read)
            new_insertion.resolve_insertion_site()
            self.resolve_insertion_site()
            return new_insertion

    def print_sum(self):
        start = self.up_closest if self.at_up_edge else f">{self.up_closest}"
        end = self.down_closest if self.at_down_edge else f"<{self.down_closest}"
        tsd = "possible_TSD" if self.up_tsd and self.down_tsd else "no_TSD"
        inconsistent = self.get_inconsistent_up_reads() + self.get_inconsistent_down_reads()
        sum_line = f"{self.chrom}\t{start}\t{end}\t{self.insertion_strand}\t{self.get_read_number()}\t{len(inconsistent)}\t{tsd}\n"
        return sum_line

    def print_reads(self, option):  # option can be "all", "consistent", "inconsistent"
        if option == "all":
            reads = self.up_reads + self.down_reads
        elif option == "consistent":
            reads = [read for read in self.up_reads + self.down_reads if not (hasattr(read, 'inconsistent') and read.inconsistent)]
        elif option == "inconsistent":
            reads = [read for read in self.up_reads + self.down_reads if hasattr(read, 'inconsistent') and read.inconsistent]
        else:
            raise ValueError(f"Unknown option: {option}")

        return "\n".join([read.print() for read in reads]) + "\n" if len(reads) > 0 else None

    def get_inconsistent_up_reads(self):
        return [read for read in self.up_reads if hasattr(read, 'inconsistent') and read.inconsistent]

    def get_inconsistent_down_reads(self):
        return [read for read in self.down_reads if hasattr(read, 'inconsistent') and read.inconsistent]

    def get_read_number(self):
        return len(self.up_reads) + len(self.down_reads)

    def within_range(self, read, distance):
        if self.chrom == read.chrom_1 and self.insertion_strand == read.insertion_strand and read.start_1 < self.last_pos + distance:
            return True
        else:
            return False


def close_insertions(insertion1, insertion2, distance):
    if insertion1.chrom == insertion2.chrom and abs(insertion1.last_pos - insertion2.last_pos) < distance:
        return True
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
            if len(insertions)> 0 and insertions[-1].within_range(read, distance):
                insertions[-1].add_read(read)
            else:
                insertion = Insertion()
                insertion.add_read(read)
                insertions.append(insertion)
    print(f"Total initial insertions: {len(insertions)}")
    new_insertions = []
    for insertion in insertions:
        insertion.resolve_insertion_site()
        if len(insertion.get_inconsistent_up_reads()) > 2:
            new_insertion = insertion.find_secondary_up_insertion()
            new_insertions.append(new_insertion)
        if len(insertion.get_inconsistent_down_reads()) > 2:
            new_insertion = insertion.find_secondary_down_insertion()
            new_insertions.append(new_insertion)

    insertions.extend(new_insertions)
    print(f"Total initial insertions: {len(insertions)}")
    with open(out, 'w') as f:
        for insertion in insertions:
            sum_line = insertion.print_sum()
            f.write(sum_line)
    
    with open(failed, 'w') as f:
        for insertion in insertions:
            if not insertion.print_reads("inconsistent") == None:
                f.write(insertion.print_reads("inconsistent"))


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