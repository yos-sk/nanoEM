#! /usr/bin/env python

'''
C → T or G → A conversion of reference genome
'''

import sys
import gzip
import os
from multiprocessing import Pool

def conversion_seq(chrom:str, sequence:str) -> tuple:

    out = dict()
    out_ct = ''
    out_ga = ''

    for base in sequence:
        if base == 'C':
            out_ct += 'T'
            out_ga += base
        elif base == 'c':
            out_ct += 't'
            out_ga += base
        elif base == 'G':
            out_ct += base
            out_ga += 'A'
        elif base == 'g':
            out_ct += base
            out_ga += 'a'
        else:
            out_ct += base
            out_ga += base

    return (chrom, out_ct, out_ga)

def conversion_ref(ref:str, out_path:str, thread:int) -> None:
    if not os.path.exists(ref):
        print("Reference file does not exit!", file=sys.stderr)
        sys.exit(1)

    if not os.path.exists(out_path):
        print("Output directory does not exit!", file=sys.stderr)
        sys.exit(1)

    ref_seq = list()

    if ref[-3:] == '.gz':
        f = gzip.open(ref, 'rt')
    elif ref[-3:] == '.fa':
        f = open(ref, 'r')
    else:
        print('Given file format is not supported', file=sys.stderr)
        sys.exit(1)

    seq = ''
    chrom = ''
    for line in f:
        if line[0] == '>':
            if chrom != '':
                ref_seq.append((chrom, seq))
            seq = ''
            chrom = line.rstrip('\n')
        else:
            seq += line
    
    ref_seq.append((chrom, seq))

    f.close()


    out_dic = dict()
    p = Pool(thread)
    results = p.starmap(conversion_seq, ref_seq)
    p.close()
    p.join()

    for result in results:
        out_dic[result[0]] = {"CT": result[1], "GA": result[2]}

    out = open(out_path + '/' + 'converted_ref.fa', 'w')

    for key in sorted(out_dic.keys()):
        print(key + "CT", file=out)
        print(out_dic[key]["CT"], end='', file=out)

    for key in sorted(out_dic.keys()):
        print(key + "GA", file=out)
        print(out_dic[key]["GA"], end = '', file=out)

    out.close()


def main():
    ref = sys.argv[1]
    out_path = sys.argv[2]
    thread = int(sys.argv[3])
    conversion_ref(ref, out_path, thread)

if __name__ == '__main__':
    main()