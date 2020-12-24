'''
C → T or G → A conversion of reference genome
'''

import sys
import gzip

def conversion(ref, flag):
    if ref[-3:] == '.gz':
        f = gzip.open(ref, 'rt')
    elif ref[-3:] == '.fa':
        f = open(ref, 'r')
    else:
        print('File format is wrong', file=sys.stderr)
        sys.exit(1)

    for line in f:
        if line[0] == '>':
            print(line.rstrip('\n') + flag)
        else:
            if flag == 'CT':
                for base in line:
                    if base == 'C':
                        print('T', end='')
                    elif base == 'c':
                        print('t', end='')
                    else:
                        print(base, end='')
            elif flag == 'GA':
                for base in line:
                    if base == 'G':
                        print('A', end='')
                    elif base == 'g':
                        print('a', end='')
                    else:
                        print(base, end='')
    f.close()

def main():
    ref = sys.argv[1]

    conversion(ref, "CT")
    conversion(ref, "GA")

if __name__ == '__main__':
    main()
