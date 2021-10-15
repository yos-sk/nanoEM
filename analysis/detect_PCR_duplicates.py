#! /usr/bin/env python

'''
Detect PCR dulicates from the 5 prime positions and the 3 prime positions of reads in a BAM file. 
'''

import argparse
import sys
import gzip
import pysam
import numpy as np
from collections import Counter

class alignment_info():
    '''
    read_name       : Read ID 
    start           : 5 prime position of an alignment (int)
    end             : 3 prime position of an alignment (int)
    strand          : Mapping orientation (str)
    is_supplementary: Whether an alignment is supplementary or not (True/False)
    '''
    def __init__(self):
        self.read_name = ''
        self.start = -1
        self.end = -1
        self.strand = ''
        self.is_supplementary = False

    def insert(self, alignment:pysam.libcalignedsegment.AlignedSegment) -> None:
        self.read_name = alignment.query_name
        self.start = alignment.reference_start
        self.end = alignment.reference_end
        self.strand = '-' if alignment.is_reverse else '+'
        self.is_supplementary = alignment.is_supplementary


def clustering_alignments(bam:str, margin:int) -> dict:
    '''
    Requirement: Sorted BAM file
    '''

    aln_info = dict()

    bamfile = pysam.AlignmentFile(bam, 'rb')
    for aln in bamfile.fetch():

        tmp_info = alignment_info()
        tmp_info.insert(aln)

        if aln.reference_name in aln_info:
            aln_info[aln.reference_name].append(tmp_info)

        else:
            aln_info[aln.reference_name] = [tmp_info]


    bamfile.close()


    out = dict()

    for chrom in aln_info:
        out[chrom] = list()
        start = 0
        end = 0
        prev_start = 0
        tmp_list = list()
        for item in aln_info[chrom]:
            if len(tmp_list) == 0:
                start = item.start
                end = item.end
                prev_start = item.start
                tmp_list = [item]
                continue

            # Confirm whether aln_info is sorted or not
            if prev_start - item.start > 0:
                print("List of alignment information is not sorted.", prev_start, item.start, sep='\t', file=sys.stderr)
                sys.exit(1)


            if abs(start - item.start) < margin and abs(end - item.end) < margin:
                #print(start, item.start, end, item.end)
                start = np.mean(np.array([i.start for i in tmp_list])) 
                end = np.mean(np.array([i.end for i in tmp_list]))
                prev_start = item.start
                tmp_list.append(item)

            else:
                if len(tmp_list) > 1:
                    out[chrom].append(tmp_list)

                start = 0
                end = 0
                prev_start = item.start
                tmp_list = list()

    return out

def validation_cluster(out_cluster:dict) -> None:
    '''
    Check any supplementary alignments
    '''

    candidates_chimera = list()
    count = 0
    count1 = 0

    for chrom in out_cluster:
        #print(out_cluster[chrom])
        for items in out_cluster[chrom]:
            supplementary_count = Counter([i.is_supplementary for i in items])

            if supplementary_count[True] > 0:
                candidates_chimera.append((chrom, items))
                continue

            print(chrom, np.mean(np.array([i.start for i in items])), np.mean(np.array([i.end for i in items])), sep='\t', end ='\t')
            count += len(items)
            if chrom != 'chrM' and chrom != 'chrEBV':
                count1 += len(items)
            for i in items:
                print(i.read_name, end=',')
            print('\tPD')

    print("# of PCR duplicates: ", count, file=sys.stderr)
    print('# of PCR duplicates (chr1-22, X)', count1, file=sys.stderr)

    if len(candidates_chimera) > 0:
        #print(candidates_chimera[0])
        check_chimera(candidates_chimera)


def check_chimera(candidates_chimera:list) -> None:
    '''
    Check PCR chimera
    '''

    chimera_count = 0
    chimera_count1 = 0
    count = 0
    count1 = 0
    for i, items in enumerate(candidates_chimera):
        chimera_info = dict()
        for k in items[1]:
            #print(k)
            chimera_info[k.read_name] = [k]

        chimera_read = set()
        for j in range(i+1, len(candidates_chimera)):
            for l in candidates_chimera[j][1]:
                #print(l)
                if l.read_name in chimera_info:
                    chimera_read.add(l.read_name)
                    chimera_info[l.read_name].append((candidates_chimera[j][0], l))

        if len(chimera_read) > 1:
            for rname in chimera_info:
                if rname not in chimera_read: continue

                pc_pos = list()

                for item in chimera_info[rname]:
                    if len(pc_pos) == 0:
                        pc_pos.append((items[0], item.start, item.end))

            flg = False
            for pos in pc_pos:
                print(pos[0], pos[1], pos[2], sep=',', end='\t')
                if pos[0] == 'chrM' or pos[0] == 'chrEBV':
                    flg = True

            for name in chimera_read:
                chimera_count += 1
                if not flg: chimera_count1 += 1
                print(name, end=',')
            print('\tPC')

        else:
            read_count = Counter([i.read_name for i in items[1]])
            se_read = list()
            for name in read_count:
                if read_count[name] > 1:
                    chimera_count += 1
                    if items[0] != 'chrM' and items[0] != 'chrEBV':
                        chimera_count1 += 1
                    se_read.append(name)
                    print(items[0], np.mean(np.array([i.start for i in items[1]])), np.mean(np.array([i.end for i in items[1]])), name, 'PC', sep='\t')

            t_flg = False 
            for i in items[1]:
                if i.read_name not in se_read:
                    if not t_flg:
                        print(items[0], np.mean(np.array([i.start for i in items[1]])), np.mean(np.array([i.end for i in items[1]])), sep='\t', end='\t')
                        print(i.read_name, end=',')
                        t_flg = True
                        count += 1
                        if items[0] != 'chrM' and items[0] != 'chrEBV':
                            count1 += 1
                        continue

                    if t_flg: print(i.read_name, end=',')
                        
                    count += 1
                    if items[0] != 'chrM' and items[0] != 'chrEBV':
                        count1 += 1
                    
            if t_flg: print('\tPCC')

    
    print('# of PCR chimeras: ', chimera_count, file=sys.stderr)
    print('# of PCR chimeras (chr1-22, X): ', chimera_count1, file=sys.stderr)
    print('# of PCR chimera candidates: ', count, file=sys.stderr)
    print('# of PCR chimera candidates (chr1-22, X): ', count1, file=sys.stderr)

                    

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("--input", "-i", required=True, help="BAM file of read")
    parser.add_argument("--M", "-m", type=int, default=10, required=True, help='')
    args = parser.parse_args()

    out = clustering_alignments(args.input, args.M)
    validation_cluster(out)




if __name__ == "__main__":
    main()


   


