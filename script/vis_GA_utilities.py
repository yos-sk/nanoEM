#! /usr/bin/env python

'''
Extract best unique alignemnt from four BAM files
Two BAM files:
1. read (C->T)
2. read (G->A)
'''

import sys
import pysam
import argparse

def complementary(seq):
    out = ''
    for base in seq[::-1]:
        if base == 'A':
            out += 'T'
        elif base == 'C':
            out += 'G'
        elif base == 'G':
            out += 'C'
        elif base == 'T':
            out += 'A'

    return out

def man_vis(bamfile):

    samfile = pysam.AlignmentFile(bamfile, 'rb')

    print(samfile.header, end='')

    for read in samfile.fetch():
        #seq = complementary(read.query_sequence)
        if read.is_reverse:
            flag = read.flag - 16
        else:
            flag = read.flag + 16

        '''
        qual = ''
        for q in read.query_qualities[::-1]:
            qual += chr(q + 33)
        '''
        q_name = read.qname
        pos = str(read.reference_start + 1)
        r_name = read.reference_name
        '''
        cigar = ''

        for c in read.cigartuples[::-1]:
            if c[0] == 0:
                cig_str = 'M'
            elif c[0] == 1:
                cig_str = 'I'
            elif c[0] == 2:
                cig_str = 'D'
            elif c[0] == 3:
                cig_str = 'N'
            elif c[0] == 4:
                cig_str = 'S'
            elif c[0] == 5:
                cig_str = 'H'
            elif c[0] == 6:
                cig_str = 'P'
            elif c[0] == 7:
                cig_str = '='
            elif c[0] == 8:
                cig_str = 'X'
            elif c[0] == 9:
                cig_str = 'B'

            cigar += str(c[1])
            cigar += cig_str
        '''
        seq = read.query_sequence
        mapq = str(read.mapping_quality)
        cigar = read.cigarstring
        NM = 'NM:i:'
        for tag in read.get_tags():
            if tag[0] == 'NM':
                NM = NM + str(NM[1])
                break
        qual = ''
        for q in read.query_qualities:
            qual += chr(q + 33)


        print('\t'.join((q_name, str(flag), r_name, pos, mapq, cigar, '*', '0', '0', seq, qual, NM)))


    samfile.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam", "-b", help="BAM file of read (G->A)")
    args = parser.parse_args()


    
    man_vis(args.bam)


if __name__ == '__main__':
    main()
