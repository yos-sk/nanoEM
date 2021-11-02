#! /usr/bin/env python

'''
Detect and remove PCR dulicates from the 5 prime positions and the 3 prime positions of reads in a BAM file.
'''

import argparse
import sys
import pysam
import numpy as np
import gzip
import subprocess

def calculate_read_length(in_fh:str) -> dict:

    out = dict()
    if in_fh[-3:] == '.gz':
        f = gzip.open(in_fh, 'rt')
    else:
        f = open(in_fh, 'r')

    rname = ''
    for i, line in enumerate(f):
        if i % 4 == 0:
            rname = line.rstrip('\n').split(' ')[0][1:]
        if i % 4 == 1:
            out[rname] = len(line.rstrip('\n'))

    f.close()

    return out

def detect_ends(in_fh:str) -> list:
    '''
    Clustering alignments with the same start and end positions in given margin.
    '''

    out = list()
    f = gzip.open(in_fh, 'rt')

    for line in f:
        items = line.rstrip('\n').split('\t')
        rname = items[0]
        if len(items[1:]) == 1:
            chrom = items[1].split(',')[0]
            start = items[1].split(',')[1]
            end = items[1].split(',')[2]
            out.append([rname, chrom + ':' + start, chrom + ':' + end])
            continue

        start_chrom = items[1].split(',')[0]
        start_start = items[1].split(',')[1]
        start_end = items[1].split(',')[2]
        end_chrom  = items[-1].split(',')[0]
        end_start = items[-1].split(',')[1]
        end_end = items[-1].split(',')[2]
        if len(items) > 2:
            is_split = True
        else:
            is_split = False

        if start_chrom < end_chrom:
            out.append([rname, start_chrom + ':' + start_start, end_chrom + ':' + end_end, is_split])
        elif end_chrom > start_chrom:
            out.append([rname, end_chrom + ':' + end_start, start_chrom + ':' + start_end, is_split])
        else:
            if start_start <=  end_start:
                out.append([rname, start_chrom + ':' + start_start, end_chrom + ':' + end_end, is_split])
            else:
                out.append([rname, end_chrom + ':' + end_start, start_chrom + ':' + start_end, is_split])

    return out


def make_cluster(ends:dict, margin:int) -> dict:
    '''
    PCR duplicates -> 1. Both of reference start and reference end are aligned
                      2. Raw read lengths are nearly same (< 250 bp).
    '''

    duplicate_cluster = dict() # key: position (chrA:xxxx-yyyy), value: list of read names

    prev_start_chrom = ''
    prev_start = 0
    tmp_list = list()
    for i, item in enumerate(sorted(ends, key = lambda x:x[1])):

        if i == 0:
            prev_start_chrom = item[1].split(':')[0]
            prev_start = int(item[1].split(':')[1])
            tmp_list = [item]
            continue

        chrom = item[1].split(':')[0]
        start = int(item[1].split(':')[1])

        if chrom == prev_start_chrom and abs(start - prev_start) < margin:
            tmp_list.append(item)
        if (i == len(ends) - 1 or chrom != prev_start_chrom or abs(start - prev_start) >= margin) and len(tmp_list) > 1:
            prev_end, end = 0, 0
            prev_end_chrom, end_chrom = '', ''
            tmp_cluster = list()
            for j, t in enumerate(sorted(tmp_list, key = lambda x:x[2])):
                if j == 0:
                    prev_end_chrom = t[2].split(':')[0]
                    prev_end = int(t[2].split(':')[1])
                    tmp_cluster = [t]
                    continue

                end_chrom = t[2].split(':')[0]
                end = int(t[2].split(':')[1])

                if end_chrom == prev_end_chrom and abs(end - prev_end) < margin:
                    tmp_cluster.append(t)
                if (j == len(tmp_list) - 1 or end_chrom != prev_end_chrom or abs(end - prev_end) >= margin) and len(tmp_cluster) > 1:
                    s_chrom = tmp_cluster[0][1].split(':')[0]
                    e_chrom = tmp_cluster[0][2].split(':')[0]
                    c_start = np.mean(np.array([int(k[1].split(':')[1]) for k in tmp_cluster]))
                    c_end  = np.mean(np.array([int(k[2].split(':')[1]) for k in tmp_cluster]))
                    duplicate_cluster[s_chrom + ':' + str(c_start) + ',' + e_chrom + ':' + str(c_end)] = [k[0] for k in tmp_cluster]
                    if j != len(tmp_list) - 1:
                        tmp_cluster = [t]
                prev_end_chrom = end_chrom
                prev_end = end

            tmp_list = [item]

        prev_start_chrom = chrom
        prev_start = start

    return duplicate_cluster

def mark_duplicates(ends:dict, duplicate_cluster:dict, in_bam:str, in_fastq:str) -> None:
    '''
    Count PCR duplicates and PCR chimera candidates
    PCR duplicates         : Reads in duplicates_cluster from make_cluster
    PCR chimera candidates : Reads in duplicates_cluster if is_split == True
    '''

    readlen_dict = calculate_read_length(in_fastq)
    mark_rnames = set()
    count = 0
    for key in duplicate_cluster:
        tmp_cluster = list()
        rname = duplicate_cluster[key][0]
        tmp_cluster.append(readlen_dict[rname])
        #print(key, rname, sep='\t', end = '')
        for i, name in enumerate(duplicate_cluster[key][1:]):
            flag = False
            for length in tmp_cluster:
                if abs(length - readlen_dict[name]) > 250:
                    flag = True
                else:
                    flag = False
                    break

            if flag:
                tmp_cluster.append(readlen_dict[name])
            else:
                mark_rnames.add(name)
                #print('\t' + name, end = '')
        #print('')
    print("# of PCR duplicates:", len(mark_rnames), sep='\t', file=sys.stderr)

    bamfile = pysam.AlignmentFile(in_bam, 'rb')
    out_bam = pysam.AlignmentFile(in_bam[:-4] + 'markdup.bam', 'wb', template = bamfile)

    for read in bamfile.fetch():
        if read.query_name in mark_rnames:
            read.flag += 1024
        out_bam.write(read)

    bamfile.close()
    subprocess.check_call(['samtools', 'index', in_bam[:-4] + 'markdup.bam'])

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--input", "-i", required=True, help="alignment information")
    parser.add_argument("--bam", "-b", required=True, help="BAM file of read")
    parser.add_argument("--M", "-m", type=int, default=50, required=False, help='')
    parser.add_argument("--fastq", "-f", required=True, help='')
    args = parser.parse_args()

    ends = detect_ends(args.input)
    duplicates = make_cluster(ends, args.M)
    mark_duplicates(ends, duplicates, args.bam, args.fastq)

if __name__ == '__main__':
    main()
