#! /usr/bin/env python

'''
C → T or G → A conversion of fastqerence genome
'''

import sys, os
import gzip
from multiprocessing import Pool
import time
from typing import Tuple

class read_fastq:
    def __init__(self):
        self.read_id = ''
        self.seq = ''
        self.plus = '+'
        self.quality = ''
        self.seq_ct = ''
        self.seq_ga = ''

    def add_read_id(self, read_id:str) -> None:
        self.read_id = read_id

    def add_seq(self, seq:str) -> None:
        self.seq = seq

    def add_quality(self, quality:str) -> None:
        self.quality = quality

def conversion_base(info:read_fastq) -> read_fastq:
    out_seq_ct = ''
    out_seq_ga = ''
    if info.seq_ct != '' and info.seq_ga != '':
        print('Duplicates in the given fastq file', file=sys.stderr)
        sys.exit(1)

    for base in info.seq:
        if base == 'C':
            out_seq_ct += 't'
        else:
            out_seq_ct += base

        if base == 'G':
            out_seq_ga += 'a'
        else:
            out_seq_ga += base

        info.seq_ct = out_seq_ct
        info.seq_ga = out_seq_ga

    return info

def split_fastq(fastq:str) -> Tuple[set, str]:
    if fastq[-3:] == '.gz':
        f = gzip.open(fastq, 'rt')
        prefix = fastq[:-6]
    elif fastq[-3:] == '.fq':
        f = open(fastq, 'r')
        prefix = fastq[:-3]
    else:
        print('File format is wrong', file=sys.stderr)
        sys.exit(1)

    num_reads = 0
    out_set = {prefix + '_0.fq'}
    w = open(prefix + '_0.fq', 'w')
    index_file = open(prefix + '_index.txt', 'w')
    for i, line in enumerate(f):
        if num_reads % 10000 == 1:
            if num_reads != 0 and i % 4 == 0:
                w.close()
                w = open(prefix + '_' + str(num_reads//10000) + '.fq', 'w')
                out_set.add(prefix + '_' + str(num_reads//10000) + '.fq')
                #print(prefix + '_' + str(num_reads//10000) + '.fq', file=sys.stderr)

        if i % 4 == 0:
            num_reads += 1
            rname = line.rstrip('\n').split(" ")[0][1:]
            print(rname, prefix + '_' + str(num_reads//10000) + '.fq', sep='\t', file=index_file)

        print(line.rstrip('\n'), file=w)

    w.close()
    index_file.close()
    assert len(out_set) == 994, '# of files conslits to the return value.'

    return out_set, prefix

def conversion_reads(fastq:str, thread:int) -> None:

    files, prefix = split_fastq(fastq)

    w1 = open(prefix + '_CT.fq', 'w')
    w2 = open(prefix + '_GA.fq', 'w')

    for k, in_fh in enumerate(files):
        seq_list = list()
        tmp_info = read_fastq()
        start = time.time()
        f = open(in_fh, 'r')
        for i, line in enumerate(f):
            if i % 4 == 0:
                if tmp_info.read_id != '' or tmp_info.seq != '' or tmp_info.quality != '':
                    print('Read info is out of format', file=sys.stderr)
                    sys.exit(1)
                tmp_info.add_read_id(line.rstrip('\n'))
            elif i % 4 == 1:
                tmp_info.add_seq(line.rstrip('\n'))
            elif i % 4 == 3:
                tmp_info.add_quality(line.rstrip('\n'))
                seq_list.append(tmp_info)

                if tmp_info.read_id == '' or tmp_info.seq == '' or tmp_info.quality == '':
                    print('Read info is out of format', file=sys.stderr)
                    sys.exit(1)

                tmp_info = read_fastq()

        f.close()
        end = time.time()
        print(k, 'Load fastq file:', end - start, sep='\t', file=sys.stderr)


        start = time.time()
        results = list()
        p = Pool(thread)
        tmp_results = p.map(conversion_base, seq_list)
        results.extend(tmp_results)
        p.close()
        p.join()
        end = time.time()

        print(k, 'Conversion of sequence:', end - start, sep='\t', file=sys.stderr)

        start = time.time()

        for result in results:
            print(result.read_id, file=w1)
            print(result.seq_ct, file=w1)
            print('+', file=w1)
            print(result.quality, file=w1)

            print(result.read_id, file=w2)
            print(result.seq_ga, file=w2)
            print('+', file=w2)
            print(result.quality, file=w2)

        end = time.time()
        print(k, 'Writing:', end - start, sep='\t', file=sys.stderr)

    w1.close()
    w2.close()


def main():
    fastq = sys.argv[1]
    thread = int(sys.argv[2])
    conversion_reads(fastq, thread)

if __name__ == '__main__':
    main()