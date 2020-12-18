'''
C → T or G → A conversion of fastqerence genome
'''

import sys
import gzip

def conversion(fastq):
    if fastq[-3:] == '.gz':
        f = gzip.open(fastq, 'rt')
    elif fastq[-3:] == '.fq':
        f = open(fastq, 'r')
    else:
        print('File format is wrong', file=sys.stderr)
        sys.exit(1)

    w1 = gzip.open('1d_pass_CT.fq.gz', 'wt')
    w2 = gzip.open('1d_pass_GA.fq.gz', 'wt')

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

def main():
    fastq = sys.argv[1]
    conversion(fastq)

if __name__ == '__main__':
    main()
