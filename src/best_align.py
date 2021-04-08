'''
Extract best unique alignemnt from four BAM files
Two BAM files:
1. read (C->T)
2. read (G->A)
'''

import sys
import gzip
import re
import pysam
import hashlib
import argparse

class cigar_string:
    def __init__(self):
        self.M = 0
        self.I = 0
        self.D = 0
        self.E = 0
        self.X = 0
        self.NM = 0

    def increment(self, cigar, count):
        if cigar == 0:
            self.M += count
        elif cigar == 1:
            self.I += count
        elif cigar == 2:
            self.D += count
        elif cigar == 7:
            self.E += count
        elif cigar == 8:
            self.X += count
        elif cigar == 9:
            self.NM += count

    def alignment_length(self):
        return self.E + self.X + self.I + self.D

    def identity_numerator(self):
        return self.E

def read_fastq(fastq):
    out = dict()
    if fastq[-3:] == '.gz':
        f = gzip.open(fastq, 'rt')
    else:
        f = open(fastq, 'r')

    for i, line in enumerate(f):
        if i % 4 == 0:
            rname = line.rstrip('\n').split()[0][1:]
            hs = int(hashlib.sha224(rname.encode()).hexdigest(), 16) % 6091
        if i % 4  == 1:
            seq = line.rstrip('\n')
            if hs in out:
                out[hs][rname] = seq
            else:
                out[hs] = dict()
                out[hs][rname] = seq

    f.close()

    return out

def extract_information(bamFile, flag, out):
    bamfile = pysam.AlignmentFile(bamFile, "rb")

    for read in bamfile.fetch():
        '''
        Filtering
        '''
        if read.flag & 4 == 4:
            continue
        if read.flag & 256 == 256:
            continue
        if read.reference_name == 'chrEBV':
            continue
        if read.mapping_quality < 20:
            continue

        if read.query_name in out:
            if flag in out[read.query_name]:
                out[read.query_name][flag].append(read)
            else:
                out[read.query_name][flag] = [read]
        else:
            out[read.query_name] = dict()
            out[read.query_name][flag] = [read]

    bamfile.close()

def complementary(seq):
    out = ''
    for base in seq[::-1]:
        if base == 'A':
            out += 'T'
        elif base == 'C':
            out += 'G'
        elif base == 'G':
            out += 'C'
        elif base == 'T':
            out += 'A'

    return out

def compare(out, bamfile, fastq):
    '''
    Compare the mapping results and decide the best unique alignment
    1: read (C->T)
    2: read (G->A)
    '''
    seq = read_fastq(fastq)
    w1 = open('output_CT.sam', 'w')
    w2 = open('output_GA.sam', 'w')

    '''
    Customize header of SAM files 
    '''
    samfile = pysam.AlignmentFile(bamfile, 'rb')
    header = str(samfile.header)
    line = header.rstrip('\n').split('\n')
    out_header = ''
    for items in line:
        if items[:3] != "@SQ":
            out_header += items
            out_header += "\n"
        elif items == "":
            continue
        else:
            t = ""
            i = items.split('\t')
            for c in i:
                if c == "@SQ":
                    t += c
                    t += '\t'
                elif c[:2] == "LN":
                    if t == "":
                        continue
                    t += c
                    t += '\n'
                elif c[-2:] == "CT":
                    t += c[:-2]
                    t += '\t'
                else:
                    t = ""
            out_header += t
    print(out_header, end='', file=w1)
    print(out_header, end='', file=w2)

    samfile.close()

    for key in out:
        '''
        Filter alignments by strand as following
        Reference    Read    Strand
        C->T         C->T      +
        G->A         C->T      -
        C->T         G->A      -
        G->A         G->A      +
        '''

        '''
        out (input dictionary of this function) format
        Key          Description
        1.           Ref. (C->T, G->A), Read (C->T)
        2.           Ref. (C->T, G->A), Read (G->A)
        '''
        aln = {1:list(), 2:list()}
        for flag in out[key]:
            for item in out[key][flag]:
                if flag == 1:
                    if not item.is_reverse and item.reference_name[-2:] == 'CT':
                        aln[flag].append(item)
                    elif item.is_reverse and item.reference_name[-2:] == 'GA':
                        aln[flag].append(item)
                else:
                    if item.is_reverse and item.reference_name[-2:] == 'CT':
                        aln[flag].append(item)
                    elif not item.is_reverse and item.reference_name[-2:] == 'GA':
                        aln[flag].append(item)
                

        '''
        Judge C-to-T conversion read or G-to-A conversion read with mapping quality
        '''

        mapq_CT = -1 if len(aln[1]) == 0 else max([i.mapping_quality for i in aln[1]])
        mapq_GA = -1 if len(aln[2]) == 0 else max([i.mapping_quality for i in aln[2]])

        if mapq_CT == mapq_GA:
            #print("This alignment is not unique", file=sys.stderr)
            continue

        f = "CT" if mapq_CT > mapq_GA else "GA"
        out_flag = 1 if f == "CT" else 2
        
        for item in out[key][out_flag]:
            s = ''
            hs = int(hashlib.sha224(key.encode()).hexdigest(), 16) % 6091
            if item.flag & 2048 != 2048 and item.is_reverse:
                s = complementary(seq[hs][key])
            elif item.flag & 2048 != 2048:
                s = seq[hs][key]
            elif item.flag & 2048 == 2048:
                if item.cigartuples[0][0] == 5: 
                    first = item.cigartuples[0][1]
                else:
                    first = 0

                if item.cigartuples[-1][0] == 5: 
                    last = item.cigartuples[-1][1]
                else:
                    last = 0

                if item.is_reverse:
                    s = complementary(seq[hs][key][last:len(seq[hs][key]) - first])
                else:
                    s = seq[hs][key][first:len(seq[hs][key]) - last]
            q_name = item.qname
            flag = str(item.flag)
            pos = str(item.reference_start + 1)
            r_name = item.reference_name[:-2]
            cigar = item.cigarstring
            mapq = str(item.mapping_quality)
            NM = 'NM:i:'
            for tag in item.get_tags():
                if tag[0] == 'NM':
                    NM = NM + str(tag[1])
                    break

            qual = ''
            for q in item.query_qualities:
                qual += chr(q + 33)

            if f == "CT": 
                print('\t'.join((q_name, flag, r_name, pos, mapq, cigar, '*', '0', '0', s, qual, NM)), file=w1)
            else:
                print('\t'.join((q_name, flag, r_name, pos, mapq, cigar, '*', '0', '0', s, qual, NM)), file=w2)


    w1.close()
    w2.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam1", "-b1", help="BAM file of read (C->T)")
    parser.add_argument("--bam2", "-b2", help="BAM file of read (G->A)")
    parser.add_argument("--fastq", "-f", help="Original fastq file")
    args = parser.parse_args()

    out = dict()

    
    extract_information(args.bam1, 1, out)
    extract_information(args.bam2, 2, out)

    compare(out, args.bam1, args.fastq)



if __name__ == '__main__':
    main()
