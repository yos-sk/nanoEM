'''
Calling methylation from the result of samtools mpileup.
'''

import sys
import argparse
import gzip

def call_methylation_forward(ifile, out):

    f = open(ifile, 'r')
    for line in f:
        items = line.rstrip('\n').split('\t')
        key = items[0] + ':' + items[1]
        base = items[2]
        num = int(items[3])
        cont = items[4]

        cnt = -1 
        flg = 0
        flag = -1
        for s in cont:
            if s == '*' or s == '#':
                continue
            if s == '^':
                flag = 0
                continue
            if flag == 0:
                flag = -1
                continue
            if s == '+' or s == '-':
                cnt = 0
                continue
            if cnt == 0:
                flg = int(s)
                cnt += 1
                continue
            if cnt > 0:
                if cnt < flg:
                    cnt += 1
                    continue
                elif cnt == flg:
                    cnt = -1
                    continue
            if s == '.' and base == 'C':
                if key in out:
                    out[key]['methyl'] += 1
                else:
                    out[key] = dict()
                    out[key]['methyl'] = 1
                    out[key]['unmethyl'] = 0
            if s == ',' and base == 'G':
                if key in out:
                    out[key]['methyl'] += 1
                else:
                    out[key] = dict()
                    out[key]['methyl'] = 1
                    out[key]['unmethyl'] = 0
            if s == 'T' and base == 'C':
                if key in out:
                    out[key]['unmethyl'] += 1
                else:
                    out[key] = dict()
                    out[key]['methyl'] = 0
                    out[key]['unmethyl'] = 1
            if s == 'a' and base == 'G':
                if key in out:
                    out[key]['unmethyl'] += 1
                else:
                    out[key] = dict()
                    out[key]['methyl'] = 0
                    out[key]['unmethyl'] = 1
    f.close()

def call_methylation_reverse(ifile, out):

    f = open(ifile, 'r')
    for line in f:
        items = line.rstrip('\n').split('\t')
        key = items[0] + ':' + items[1]
        base = items[2]
        num = int(items[3])
        cont = items[4]

        cnt = -1
        flg = 0
        flag = -1
        for s in cont:
            if s == '*' or s == '#':
                continue
            if s == '^':
                flag = 0
                continue
            if flag == 0:
                flag = -1
                continue
            if s == '+' or s == '-':
                cnt = 0
                continue
            if cnt == 0:
                flg = int(s)
                cnt += 1
                continue
            if cnt > 0:
                if cnt < flg:
                    cnt += 1
                    continue
                elif cnt == flg:
                    cnt = -1
                    continue
            if s == '.' and base == 'G':
                if key in out:
                    out[key]['methyl'] += 1
                else:
                    out[key] = dict()
                    out[key]['methyl'] = 1
                    out[key]['unmethyl'] = 0
            if s == ',' and base == 'C':
                if key in out:
                    out[key]['methyl'] += 1
                else:
                    out[key] = dict()
                    out[key]['methyl'] = 1
                    out[key]['unmethyl'] = 0
            if s == 'A' and base == 'G':
                if key in out:
                    out[key]['unmethyl'] += 1
                else:
                    out[key] = dict()
                    out[key]['methyl'] = 0
                    out[key]['unmethyl'] = 1
            if s == 't' and base == 'C':
                if key in out:
                    out[key]['unmethyl'] += 1
                else:
                    out[key] = dict()
                    out[key]['methyl'] = 0
                    out[key]['unmethyl'] = 1
    f.close()

def main():
    fw = sys.argv[1]
    re = sys.argv[2]

    out = dict()

    call_methylation_forward(fw, out)
    call_methylation_reverse(re, out)

    print('chromosome', 'position', 'methyl', 'unmethyl', 'freqeucy', sep='\t') 
    for key in out:
        c, p = key.split(':')
        freq = out[key]['methyl'] / (out[key]['methyl'] + out[key]['unmethyl']) 
        print(c, p, out[key]['methyl'], out[key]['unmethyl'], freq, sep='\t') 

if __name__ == '__main__':
    main()
