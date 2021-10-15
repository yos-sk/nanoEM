#! /usr/bin/env python

import sys

def extract_CpG(reference):
    out = list()

    f = open(reference, 'r')

    position = 0
    for line in f:
        if line[0] == '>':
            chrom = line.rstrip('\n')[1:]
            position = 0
        else:
            sequence = line.rstrip('\n')
            for i, base in enumerate(sequence):
                if i == 0:
                    continue
                #if base.upper() == 'G' and sequence[i-1].upper() == 'C':
                if base == 'G' and sequence[i-1] == 'C':
                    out.append(chrom + ':' + str(position + i))

            position += len(sequence)

    return out


def output(out):
    '''
    Output extracted positions by bed format
    '''

    for item in out:
        items = item.split(':')
        print(items[0], int(items[1]) - 1, items[1], sep='\t')
        #print(items[0], items[1], int(items[1]) + 1, sep='\t')


def main():
    out = extract_CpG(sys.argv[1])
    output(out)

if __name__ == '__main__':
    main()