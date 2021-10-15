#! /usr/bin/env python

'''
C → T or G → A conversion of fastqerence genome
'''

import sys
import gzip

def conversion_reads(fastq):
    if fastq[-3:] == '.gz':
        f = gzip.open(fastq, 'rt')
        prefix = fastq[:-6]
    elif fastq[-3:] == '.fq':
        f = open(fastq, 'r')
        prefix = fastq[:-3]
    else:
        print('File format is wrong', file=sys.stderr)
        sys.exit(1)

    w1 = gzip.open(prefix + '_CT.fq.gz', 'wt')
    w2 = gzip.open(prefix + '_GA.fq.gz', 'wt')

    for i, line in enumerate(f):
        if i % 4 != 1:
            print(line, end='', file=w1)
            print(line, end='', file=w2)
        else:
            '''
            C-to-T conversion
            '''
            for base in line:
                if base == 'C':
                    print('t', end='', file=w1)
                else:
                    print(base, end='', file=w1)
            '''
            G-to-A conversion
            '''
            for base in line:
                if base == 'G':
                    print('a', end='', file=w2)
                else:
                    print(base, end='', file=w2)

    w1.close()
    w2.close()

