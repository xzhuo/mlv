import os
import argparse
import re
from statistics import median
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
        # self.strand_2 is the virus strand, not the mouse strand.
        if self.strand_1 == '+':
            self.position = "up"
            self.edge = self.end_1
        if self.strand_1 == '-':
            self.position = "down"
            self.edge = self.start_1
        # did not filter out 5' soft clip, since we don't have the mapping strand on mouse in the bedpe file.
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
            #     self.up_edgies_softclip.append(read.edge)
            # self.up_softclip = self.up_softclip and read.soft_clip_1
            # self.up_limit = max(self.up_limit, read.edge)
            # if read.soft_clip_1:
            #     self.up_edgies.append(read.edge)
            if read.any_soft_clip:
                self.up_tsd = True

        if read.position == "down":
            self.down_reads.append(read)
            # self.down_softclip = self.down_softclip and read.soft_clip_1
            # if self.down_limit == 0:
            #     self.down_limit = read.edge
            # else:
            #     self.down_limit = min(self.down_limit, read.edge)
            # if read.soft_clip_1:
            #     self.down_edgies.append(read.edge)
            if read.any_soft_clip:
                self.down_tsd = True

    def resolve_insertion_site(self, tsd = 4, wiggle_distance = 20):
        for read in self.up_reads:
            conflict_up_num = len([i for i in self.up_reads if i.edge > read.edge + wiggle_distance])
            conflict_down_num = len([i for i in self.down_reads if i.edge < read.edge - tsd - wiggle_distance])
            conflict_num = conflict_up_num + conflict_down_num
            read.edge_score = 100 / (conflict_num + 1) # Convert the number of other reads conflict with it being at the edge to a score. 100 / (conflict_num + 1)

        for read in self.down_reads:
            conflict_up_num = len([i for i in self.up_reads if i.edge > read.edge + tsd + wiggle_distance])
            conflict_down_num = len([i for i in self.down_reads if i.edge < read.edge - wiggle_distance])
            conflict_num = conflict_up_num + conflict_down_num
            read.edge_score = 100 / (conflict_num + 1) # Convert the number of other reads conflict with it being at the edge to a score. 100 / (conflict_num + 1)
        if len(self.up_reads + self.down_reads) > 0:
            self.max_score = max([read.edge_score for read in self.up_reads + self.down_reads])
        # if len(self.down_reads) > 0:
        #     self.max_down_score = max([read.edge_score for read in self.down_reads])

        # label closest, at_edge on all reads,
        # label upstream and downstream closest position and whether the edge is soft clipped on the insertion.
        up_closest_all = []
        down_closest_all = []
        for read in self.up_reads:
            if read.edge_score == self.max_score:
                read.closest = True
                up_closest_all.append(read.edge)
                if read.soft_clip_1:
                    read.at_edge = True
                    self.at_up_edge = True
        for read in self.down_reads:
            if read.edge_score == self.max_score:
                read.closest = True
                down_closest_all.append(read.edge)
                if read.soft_clip_1:
                    read.at_edge = True
                    self.at_down_edge = True
        self.up_closest = max(up_closest_all) if len(up_closest_all) > 0 else 0
        self.down_closest = min(down_closest_all) if len(down_closest_all) > 0 else self.last_pos
        self.label_inconsistent_reads(up_pos=self.up_closest, down_pos=self.down_closest, tsd=tsd, wiggle_distance=wiggle_distance)

    def label_inconsistent_reads(self, up_pos, down_pos, tsd=4, wiggle_distance=20):
        self.reset_inconsistent()
        for read in self.up_reads:
            if up_pos and read.edge - wiggle_distance > up_pos:
                read.inconsistent = True
            if down_pos and read.edge - wiggle_distance > down_pos + tsd:
                read.inconsistent = True
        for read in self.down_reads:
            if up_pos and read.edge + wiggle_distance < up_pos - tsd:
                read.inconsistent = True
            if down_pos and read.edge + wiggle_distance < down_pos:
                read.inconsistent = True

    def reset_inconsistent(self):
        for read in self.up_reads:
            read.inconsistent = False
        for read in self.down_reads:
            read.inconsistent = False

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

    def diagnose(self, tsd = 4, wiggle_distance = 20):
        up_softclip_pos = [read.edge for read in self.up_reads if read.soft_clip_1]
        merged_up_edges = merge_pos(up_softclip_pos, wiggle_distance)
        good_up_edges = [i for i in merged_up_edges if i["count"] > 1]
        down_softclip_pos = [read.edge for read in self.down_reads if read.soft_clip_1]
        merged_down_edges = merge_pos(down_softclip_pos, wiggle_distance)
        good_down_edges = [i for i in merged_down_edges if i["count"] > 1]
        if len(good_up_edges) == 0 or len(good_down_edges) == 0:
            pass
        elif len(good_up_edges) == 1 and len(good_down_edges) == 1:
            up_edge = good_up_edges[0]['median_pos']
            down_edge = good_down_edges[0]['median_pos']
            tsd_distance = abs(down_edge - up_edge)
            print(f"Found insertion at {self.chrom}:{up_edge}-{down_edge} with high TSD {tsd_distance}")
            self.reset_scores()
            self.resolve_insertion_site(tsd=tsd_distance)
        else:
            if len(good_up_edges) > 1 or len(good_down_edges) > 1:
                # since the wiggle_distance is far larger than the tsd, directly compare upstream and downstream edges by merging:
                merged_combined_edges = merge_pos([u["median_pos"] for u in good_up_edges] + [d["median_pos"] for d in good_down_edges], wiggle_distance = 20)
                # Now we can find good edges that are supported by both up and down reads
                good_edges_up_down = [c for c in merged_combined_edges if c["count"] > 1]
                if len(good_edges_up_down) > 0:
                    sorted_combined_good_edges = sorted(good_edges_up_down, key=lambda i: i['count'], reverse=True)
                    self.reset_scores()
                    best_edge = sorted_combined_good_edges[0]
                    # upstream edge should be the second position.
                    up_pos = best_edge['pos_list'][1] if best_edge['pos_list'] else None
                    down_pos = best_edge['pos_list'][0] if best_edge['pos_list'] else None
                    self.label_inconsistent_reads(up_pos=up_pos, down_pos=down_pos)
                    new_insertions = []
                    if len(self.get_inconsistent_up_reads()) > 2:
                        new_insertions.append(self.find_secondary_up_insertion())
                    if len(self.get_inconsistent_down_reads()) > 2:
                        new_insertions.append(self.find_secondary_down_insertion())
                    return new_insertions

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
            separate_reads = [read for read in self.down_reads if read.edge > new_insertion.up_closest - tsd - wiggle_distance]
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
            separate_reads = [read for read in self.up_reads if read.edge < new_insertion.down_closest + tsd + wiggle_distance]
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

    def print_reads(self, option = "all"):  # option can be "all", "consistent", "inconsistent"
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


def merge_pos(list, wiggle_distance): # input is a list of integers.
    list.sort()
    merged = [] # a list of dicts of {"pos_list": [pos], "count": count, "median_pos": median_pos}
    for item in list:
        if not merged:
            merged.append({"pos_list": [item], "count": 1})
        else:
            if abs(merged[-1]["pos_list"][-1] - item) < wiggle_distance:
                merged[-1]["pos_list"].append(item)
                merged[-1]["count"] += 1
            else:
                merged.append({"pos_list": [item], "count": 1})
    for i in merged:
        i["median_pos"] = int(median(i["pos_list"]))
    return merged

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

            new_insertion_list = insertion.diagnose()
            if new_insertion_list:
                new_insertions.extend(new_insertion_list)
        if len(insertion.get_inconsistent_down_reads()) > 2:

            new_insertion_list = insertion.diagnose()
            if new_insertion_list:
                new_insertions.extend(new_insertion_list)

    insertions.extend(new_insertions)
    insertions.sort(key=lambda x: (x.chrom, x.last_pos))
    print(f"Total initial insertions: {len(insertions)}")
    with open(out, 'w') as f:
        for insertion in insertions:
            if not (insertion.up_tsd and insertion.down_tsd):
                sum_line = insertion.print_sum()
                f.write(sum_line)

    with open(tsd, 'w') as f:
        for insertion in insertions:
            if insertion.up_tsd and insertion.down_tsd:
                sum_line = insertion.print_sum()
                bedpe = insertion.print_reads("all")
                f.write(f"# {sum_line}{bedpe}\n")


def main():
    parser = argparse.ArgumentParser(description='parse a bedpe file and summerize insertion strands')
    parser.add_argument('-i', '--input', type=str, required=True,
                        help='input bedpe file')
    parser.add_argument('-o', '--out', type=str, required=True,
                        help='output summery of no TSD records')
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