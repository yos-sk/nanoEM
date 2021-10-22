#! /usr/bin/env python

'''
Calling methylation from the result of samtools mpileup.
'''


import sys
import argparse
import gzip

def parse_pileup(pileup:str) -> list:
    base = list()
    cnt = 0
    while cnt < len(pileup):
            if pileup[cnt] in ['*' or '#']: 
                base.append(pileup[cnt])
                cnt += 1
            elif pileup[cnt] == '^': cnt += 2
            elif pileup[cnt] in ['+', '-']:
                cnt += 1
                n = ''
                while pileup[cnt] in ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']:
                    n += pileup[cnt]
                    cnt += 1
                cnt += int(n)
            elif pileup[cnt] in ['>', '<', '$']: cnt += 1
            else:
                base.append(pileup[cnt])
                cnt += 1

    return base

def base_identifier(cpgs:str) -> dict:

    out = dict()

    if cpgs[-3:] == '.gz':
        f = gzip.open(cpgs, 'rt')
    else:
        f = open(cpgs, 'r')
    
    for i, line in enumerate(f):
        items = line.rstrip('\n').split('\t')
        key = items[0] + ':' + items[2]
        if i % 2 == 0:
            out[key] = 'C'
        else:
            out[key] = 'G'
    
    f.close()
    return out

def call_methylation(in_fh:str, out:dict, strand:str, flag_ref:bool, cpgs:str):

    f = open(in_fh, 'r')
    if not flag_ref:
        base_dic = base_identifier(cpgs)
    
    for line in f:
        items = line.rstrip('\n').split('\t')
        key = items[0] + ':' + items[1]
        if flag_ref:
            base = items[2]
        else:
            base = base_dic[key]
             
        num = int(items[3])
        if num == 0: continue
        p_base = parse_pileup(items[4])

        assert len(p_base) == num,  'Number of bases: {0} conflicts with parsed pileup; {1}'.format(num, len(p_base))

        if key not in out:
            out[key] = {'methyl':0, 'unmethyl':0}
        
        if flag_ref:
            if strand == '+':
                for s in p_base:
                    if s == '.' and base == 'C':
                        out[key]['methyl'] += 1
                    if s == ',' and base == 'G':
                        out[key]['methyl'] += 1
                    if s == 'T' and base == 'C':
                        out[key]['unmethyl'] += 1
                    if s == 'a' and base == 'G':
                        out[key]['unmethyl'] += 1
            else:
                for s in p_base:
                    if s == '.' and base == 'G':
                        out[key]['methyl'] += 1
                    if s == ',' and base == 'C':
                        out[key]['methyl'] += 1
                    if s == 'A' and base == 'G':
                        out[key]['unmethyl'] += 1
                    if s == 't' and base == 'C':
                        out[key]['unmethyl'] += 1
        else:
            if strand == '+':
                for s in p_base:
                    if s == 'C' and base == 'C':
                        out[key]['methyl'] += 1
                    if s == 'g' and base == 'G':
                        out[key]['methyl'] += 1
                    if s == 'T' and base == 'C':
                        out[key]['unmethyl'] += 1
                    if s == 'a' and base == 'G':
                        out[key]['unmethyl'] += 1
            else:
                for s in p_base:
                    if s == 'G' and base == 'G':
                        out[key]['methyl'] += 1
                    if s == 'c' and base == 'C':
                        out[key]['methyl'] += 1
                    if s == 'A' and base == 'G':
                        out[key]['unmethyl'] += 1
                    if s == 't' and base == 'C':
                        out[key]['unmethyl'] += 1

    f.close()