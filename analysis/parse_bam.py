#! /usr/bin/env python

'''
Parse BAM file, and construct mapping information by each read.
'''

import pysam
import gzip
import argparse
import sys


def parse_bam(bamFile:str, mapq:int, length:int) -> None:
    ex_info = dict()

    bamfile = pysam.AlignmentFile(bamFile, 'rb')


    for read in bamfile.fetch():
        if read.flag & 4 == 4 or read.mapping_quality < mapq: continue

        if read.query_name not in ex_info:
            ex_info[read.query_name] = list()

        read_length = 0
        read_start = 1
        strand = '-' if read.is_reverse else '+'
        reference_start = read.get_reference_positions()[0] + 1
        reference_end = read.get_reference_positions()[-1] + 1
        tmp_start = reference_start
        r_len = read.infer_query_length()

        if read.cigartuples[0][0] == 4:
            r_len = r_len - read.cigartuples[0][1]

        if read.cigartuples[-1][0] == 4:
            r_len = r_len - read.cigartuples[-1][1]

        if not read.is_reverse:
            if read.cigartuples[0][0] == 4 or read.cigartuples[0][0] == 5:
                read_start += read.cigartuples[0][1]
        else:
            if read.cigartuples[-1][0] == 4 or read.cigartuples[-1][0] == 5:
                read_start += read.cigartuples[-1][1]

        read_end = read_start + r_len - 1

        '''
        read.cigartuples: (cigar, length)
        Cigar (example):
            0 -> Match
            1 -> Insertion
            2 -> Deletion
            3 -> Skipped region from the reference
            4 -> Soft clipping
            5 -> Hard clipping
        '''

        for cigar in read.cigartuples:
            if cigar[0] == 0:
                read_length += cigar[1]
                tmp_start += cigar[1]
            if cigar[0] == 1:
                read_length += cigar[1]
            if cigar[0] == 2:
                if cigar[1] > length:
                    if read.is_reverse:
                        ex_info[read.query_name].append([read.reference_name, reference_start, tmp_start - 1, strand, 
                                                         r_len - read_length + 1, read_end, read.mapping_quality])
                        read_end = r_len - read_length
                        reference_start = tmp_start + cigar[1] - 1
                    else:
                        ex_info[read.query_name].append([read.reference_name, reference_start, tmp_start - 1, strand, 
                                                         read_start, read_start + read_length-1, read.mapping_quality])
                        read_start += read_length
                        reference_start = tmp_start + cigar[1] - 1
                tmp_start += cigar[1]
            if cigar[0] == 3:
                tmp_start += cigar[1]



        if read_length != r_len:
            print('Read length is inconsistent', read.query_name, read_length, r_len, sep='\t', file=sys.stderr)
            sys.exit(1)

        ex_info[read.query_name].append([read.reference_name, reference_start, reference_end, strand, read_start, read_end, read.mapping_quality])

    bamfile.close()

    w = gzip.open('alignment_information.txt.gz', 'wt')

    for k in ex_info:
        print(k, file=w, end='\t')
        for i, item in enumerate(sorted(ex_info[k], key=lambda x:x[4])): # sorted list by read_start
            if i+1 < len(ex_info[k]):
                print(','.join([str(i) for i in item]), file=w, end='\t')
            else:
                print(','.join([str(i) for i in item]), file=w)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bam', '-b' , help="PATH of BAM file")
    parser.add_argument('--mapq', '-q', type=int, help="Threshold of mapping quality")
    parser.add_argument('--length', '-l', type=int, help="Threshold of SV length")
    args = parser.parse_args()

    parse_bam(args.bam, args.mapq, args.length)

if __name__ == '__main__':
    main()