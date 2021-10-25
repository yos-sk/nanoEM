#! /usr/bin/env python

'''
Detect and remove PCR dulicates from the 5 prime positions and the 3 prime positions of reads in a BAM file. 
'''

import argparse
import sys
import pysam
import numpy as np
from collections import Counter

class alignment_info():
    '''
    read_name       : Read ID 
    start           : 5 prime position of an alignment (int)
    end             : 3 prime position of an alignment (int)
    mapq            : Mapping quality (int)
    strand          : Mapping orientation (str)
    is_supplementary: Whether an alignment is supplementary or not (True/False)
    '''
    def __init__(self):
        self.read_name = ''
        self.start = -1
        self.end = -1
        self.mapq = -1
        self.strand = ''
        self.is_supplementary = False

    def insert(self, alignment:pysam.libcalignedsegment.AlignedSegment) -> None:
        self.read_name = alignment.query_name
        self.start = alignment.reference_start
        self.end = alignment.reference_end
        self.mapq = alignment.mapping_quality
        self.strand = '-' if alignment.is_reverse else '+'
        self.is_supplementary = alignment.is_supplementary


def clustering_alignments(bam:str, margin:int) -> dict:
    '''
    Clustering alignments with the same start and end positions in given margin.
    Requirement: Sorted BAM file and the margin
    '''

    aln_info = dict() # key: chromosome, value: a list of alignment_info classes

    bamfile = pysam.AlignmentFile(bam, 'rb')
    for aln in bamfile.fetch():

        tmp_info = alignment_info()
        tmp_info.insert(aln)

        if aln.reference_name in aln_info:
            aln_info[aln.reference_name].append(tmp_info)

        else:
            aln_info[aln.reference_name] = [tmp_info]


    bamfile.close()


    out_cluster = dict() # key: position (chrA:xxxx-yyyy), value: clustered aignment info
    out_reads = dict() # key: read name, value: clustered reads

    for chrom in aln_info:
        start, end, prev_start = 0, 0, 0
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
                start = np.mean(np.array([i.start for i in tmp_list])) 
                end = np.mean(np.array([i.end for i in tmp_list]))
                prev_start = item.start
                tmp_list.append(item)

            else:
                key = chrom + ':' + str(start) + '-' + str(end)
                if len(tmp_list) > 1:
                    out_cluster[key] = tmp_list
                    for item in tmp_list:
                        if item.read_name in out_reads:
                            out_reads[item.read_name].append(key)
                        else:
                            out_reads[item.read_name] = [key]

                start = 0
                end = 0
                prev_start = item.start
                tmp_list = list()

    return out_cluster, out_reads

def count_duplicates(out_cluster:dict, out_reads:dict) -> list():
    '''
    Count PCR duplicates, PCR chimeras, and PCR chimera candidates
    PCR duplicates         : Reads in out_cluster from clustering_alignments
    PCR chimeras           : Reads in out_reads if length of out_reads > 1
    PCR chimera candidates : Reads not in out_reads but is_supplementary is true
    '''

    count = 0
    cc_names = [rname for rname in out_reads if len(out_reads) > 1]
    cc_candidates = list()
    out = list()

    for items in out_cluster.values():
        #chrom = key.split(':')[0]
        #print(chrom, np.mean(np.array([i.start for i in items])), np.mean(np.array([i.end for i in items])), sep='\t', end ='\t')
        count += len(items)

        max_mapq = 0
        tmp_rname = ''
        for item in items:
            if item.is_supplementary and item.read_name not in cc_names: cc_candidates.append(item.read_name)
            if max_mapq < item.mapq:
                tmp_rname = item.read_name
        out.append(tmp_rname)
            

    print("# of PCR duplicates:", count, sep='\t', file=sys.stderr)
    print('# of PCR chimeras:', len(cc_names), sep='\t', file=sys.stderr)
    print('# of PCR chimera candidates:', len(cc_candidates), sep='\t', file=sys.stderr)

    return out_reads.keys(), out

def mark_dup(rnames:list, ex_rnames:list, in_bam:str) -> None:
    '''
    Mark PCR duplicates (flag: 1024) to an input bam file
    '''
    bamfile = pysam.AlignmentFile(in_bam, 'rb')
    out_bam = pysam.AlignemtFile(in_bam[:-4] + 'markdup.bam', 'wb', template = bamfile)


    for read in bamfile.fetch():
        if read.query_name in rnames and read.query_name not in ex_rnames:
            read.flag += 1024
        out_bam.write(read)

    bamfile.close()

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--input", "-i", required=True, help="BAM file of read")
    parser.add_argument("--M", "-m", type=int, default=10, required=False, help='')
    args = parser.parse_args()

    out_cluster, out_reads = clustering_alignments(args.input, args.M)
    rname_list, ex_rnames = count_duplicates(out_cluster, out_reads)
    mark_dup(rname_list, ex_rnames, args.input)

if __name__ == '__main__':
    main()
