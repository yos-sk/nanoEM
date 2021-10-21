#! /usr/bin/env python

'''
C → T or G → A conversion of reference genome
'''

import sys
import gzip
import os

def conversion_ref(ref:str, out_path:str) -> None:
    
    if not os.path.exists(ref):
        print("Reference file does not exit!", file=sys.stderr)
        sys.exit(1)
    
    if not os.path.exists(out_path):
        print("Output directory does not exit!", file=sys.stderr)
        sys.exit(1)
        
    out = gzip.open(out_path + '/' + 'converted_ref.fa.gz', 'wt')

    for flag in ["CT", "GA"]:
        if ref[-3:] == '.gz':
            f = gzip.open(ref, 'rt')
        elif ref[-3:] == '.fa':
            f = open(ref, 'r')
        else:
            print('File format is wrong', file=sys.stderr)
            sys.exit(1)

        for line in f:
            if line[0] == '>':
                print(line.rstrip('\n') + flag, file=out)
            else:
                if flag == 'CT':
                    for base in line:
                        if base == 'C':
                            print('T', end='', file=out)
                        elif base == 'c':
                            print('t', end='', file=out)
                        else:
                            print(base, end='', file=out)
                elif flag == 'GA':
                    for base in line:
                        if base == 'G':
                            print('A', end='', file=out)
                        elif base == 'g':
                            print('a', end='', file=out)
                        else:
                            print(base, end='', file=out)
        f.close()
    out.close()
