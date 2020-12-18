'''
Compare the results of methylation calling from two platforms (ex) WGBS vs nanopolish).

'''

import math
import argparse
import sys
import csv
import gzip

def make_key(c, s, e):
    return c + ":" + str(s) + "-" + str(e)

def parse_key(key):
    return key.replace("-", ":").split(":")

class MethylationStats:
    def __init__(self, num_reads, num_methylated):
        self.num_reads = num_reads
        self.num_methylated_reads = num_methylated

    def methylation_frequency(self):
        return float(self.num_methylated_reads) / self.num_reads

    def mvalue(self):
        alpha = 1
        return math.log2((self.num_methylated_reads + alpha) / (self.num_reads - self.num_methylated_reads + alpha))

def load_nanopolish(filename):
    out = dict()
    csv_reader = csv.DictReader(open(filename), delimiter='\t')

    chromosomes = ['chr1', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr2', 'chr20', 'chr21', 'chr22', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chrX', 'chrY']
    for record in csv_reader:
        if record["chromosome"] not in chromosomes:
            continue

        num_reads = int(record["called_sites"])
        methylated_reads = int(record["called_sites_methylated"])
        key = make_key(record["chromosome"], int(record["start"])+1, int(record["end"])+1)
        out[key] = MethylationStats(num_reads, methylated_reads)

        # skip non-singleton, for now
        #if int(record["num_motifs_in_group"]) > 1:
        #    continue

    return out

def load_bisulfite(filename):
    out = dict()
    chromosomes = ['chr1', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr2', 'chr20', 'chr21', 'chr22', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chrX', 'chrY']
    if filename[-3:] == '.gz':
        f = gzip.open(filename, 'rt')
    else:
        f = open(filename, 'r')
    for line in f:
        items = line.rstrip('\n').split('\t')
        chrom = items[0]
        site = int(items[2])
        methyl = int(items[4])
        unmethyl = int(items[5])
        if chrom not in chromosomes:
            continue
        num_reads = methyl + unmethyl
        key = make_key(chrom, site, site)
        out[key] = MethylationStats(num_reads, methyl)

    f.close()

    return out

def load_tsv(filename):
    out = dict()
    chromosomes = ['chr1', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr2', 'chr20', 'chr21', 'chr22', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chrX', 'chrY']
    if filename[-3:] == '.gz':
        f = gzip.open(filename, 'rt')
    else:
        f = open(filename, 'r')

    next(f)
    for line in f:
        items = line.rstrip('\n').split('\t')
        chrom = items[0]
        site = int(items[1])
        methyl = int(items[2])
        unmethyl = int(items[3])
        if chrom not in chromosomes:
            continue
        num_reads = methyl + unmethyl
        key = make_key(chrom, site, site)
        out[key] = MethylationStats(num_reads, methyl)

    f.close()

    return out


def compare_methylation(ifile1, ifile2, flag1, flag2):
    if flag1 == 'WGBS' or flag1 == 'EM':
        set1 = load_bisulfite(ifile1)
    elif flag1 == 'nanopolish':
        set1 = load_nanopolish(ifile1)
    else:
        set1 = load_tsv(ifile1)
    print('Loading of file1 data has completed', file=sys.stderr)
    if flag2 == 'WGBS' or flag2 == 'EM':
        set2 = load_bisulfite(ifile2)
    elif flag2 == 'nanopolish':
        set2 = load_nanopolish(ifile2)
    else:
        set2 = load_tsv(ifile2)
    print('Loading of file2 data has completed', file=sys.stderr)
    print("key\tdepth_1\tfrequency_1\tMvalue_1\tdepth_2\tfrequency_2\tMvalue_2")
    for key in set1:
        if key in set2:

            d1 = set1[key]
            d2 = set2[key]

            if d1.num_reads == 0 or d2.num_reads == 0:
                continue
            print(key, d1.num_reads, d1.methylation_frequency(), d1.mvalue(), d2.num_reads, d2.methylation_frequency(), d2.mvalue(), sep='\t')

        else:
            d1 = set1[key]
            print(key, d1.num_reads, d1.methylation_frequency(), d1.mvalue(), None, None, None, sep='\t')



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--file1", "-f1", help="PATH of methylation calling file1")
    parser.add_argument("--file2", "-f2", help="PATH of methylation calling file2")
    parser.add_argument("--flag1", "-l1", help="Category of file1")
    parser.add_argument("--flag2", "-l2", help="Category of file2")
    args = parser.parse_args()

    compare_methylation(args.file1, args.file2, args.flag1, args.flag2)
    
