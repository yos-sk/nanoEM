'''
Calling methylation from the result of samtools mpileup.
'''

import sys
import argparse
import gzip

def parse_pileup(pileup:str) -> list:
    base = list()
    while cnt < len(pileup):
            if pileup[cnt] in ['*' or '#']: 
                base.append(pileup[cnt])
                cnt += 1
            elif pileup[cnt] == '^': cnt += 2
            elif pileup[cnt] in ['+', '-']:
                sub_cnt = cnt + 1
                n = ''
                while pileup[sub_cnt] in ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']:
                    n += pileup[sub_cnt]
                    sub_cnt += 1
                cnt = sub_cnt + int(n)
            elif pileup[cnt] in ['>', '<', '$']: cnt += 1
            else:
                base.append(pileup[cnt])
                cnt += 1

    return base

def call_methylation(in_fh:str, out:dict, strand:str, flag_ref:bool):

    f = open(in_fh, 'r')
    prev_key = ''
    for line in f:
        items = line.rstrip('\n').split('\t')
        key = items[0] + ':' + items[1]
        if flag_ref:
            base = items[2]
        else:
            if prev_key == '':
                base = 'C'
            elif prev_key.split(':')[0] != key.split(':')[0]:
                base = 'C'
            elif int(key.split(':')[1]) - int(prev_key.split(':')[1]) == 1:
                base == 'C'
            else:
                base == 'G' 
        num = int(items[3])
        p_base = parse_pileup(items[4])

        assert len(p_base) == num,  'Number of bases: [{0}] conflicts with parsed pileup; [{1}].format(num, len(p_base))'

        if key not in out:
            out[key] = {'methyl':0, 'unmethyl':0}
        
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

    f.close()


def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("--Reference", "-r", type=bool, default=True, required=True, help='')
    args = parser.parse_args()

    fw = sys.argv[1]
    re = sys.argv[2]
    
    out = dict()
    
    call_methylation(fw, out, '+', args.Reference)
    call_methylation(re, out, '-', args.Reference)

    print('chromosome', 'position', 'methyl', 'unmethyl', 'freqeucy', sep='\t') 
    for key in out:
        c, p = key.split(':')
        freq = out[key]['methyl'] / (out[key]['methyl'] + out[key]['unmethyl']) 
        print(c, p, out[key]['methyl'], out[key]['unmethyl'], freq, sep='\t') 

if __name__ == '__main__':
    main()
