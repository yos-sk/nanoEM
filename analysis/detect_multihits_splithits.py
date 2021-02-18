'''
Detect reads with multiple hits and split hits from BAM file
'''

import sys
import argparse
import pysam
import hashlib


def detect_from_wgs(bam:str, rname:bool) -> None:

    multi_reads = set()
    split_reads = set()

    bamfile = pysam.AlignmentFile(bam, 'rb')
    
    for read in bamfile.fetch():
        if read.flag & 256 != 256 and read.flag & 2048 != 2048:
            continue
       
        if read.flag & 256 == 256:
            multi_reads.add(read.query_name)
      
        if read.flag & 2048 == 2048:
            split_reads.add(read.query_name)
   
    bamfile.close()
 
    print('# of reads with multiple hits: ', len(multi_reads) ,sep='\t', file=sys.stderr)
    print('# of reads with split hits: ', len(split_reads) ,sep='\t', file=sys.stderr)

    if rname:
        for name in multi_reads:
            print(name, 'Multiple hits', sep='\t')

        for name in split_reads:
            print(name, 'Split hits', sep='\t')

    detect_unique_alignment_wgs(multi_reads, bam)


def detect_from_nanoEM(bam_CT:str, bam_GA:str, bam_original_CT:str, bam_original_GA:str, rname:bool) -> None:

    dict_reads = {"CT": dict(), "GA": dict()}

    bamfile_CT = pysam.AlignmentFile(bam_CT, 'rb')

    for read in bamfile_CT.fetch():
        hs = int(hashlib.sha224(read.query_name.encode()).hexdigest(), 16) % 6091

        if hs in dict_reads["CT"]:
            dict_reads["CT"][hs].add(read.query_name)
        else:
            dict_reads["CT"][hs] = {read.query_name}

    bamfile_CT.close()

    bamfile_GA = pysam.AlignmentFile(bam_GA, 'rb')

    for read in bamfile_GA.fetch():
        hs = int(hashlib.sha224(read.query_name.encode()).hexdigest(), 16) % 6091

        if hs in dict_reads["GA"]:
            dict_reads["GA"][hs].add(read.query_name)
        else:
            dict_reads["GA"][hs] = {read.query_name}

    bamfile_GA.close()


    multi_reads = {"CT": set(), "GA": set()}
    split_reads = {"CT": set(), "GA": set()}
    bamfile_original_CT = pysam.AlignmentFile(bam_original_CT, 'rb')

    for read in bamfile_original_CT.fetch():
        hs = int(hashlib.sha224(read.query_name.encode()).hexdigest(), 16) % 6091

        if read.query_name not in dict_reads["CT"][hs]: continue

        if read.flag & 256 != 256 and read.flag & 2048 != 2048: continue

        if read.flag & 256 == 256:
            multi_reads["CT"].add(read.query_name)

        if read.flag & 2048 == 2048:
            split_reads["CT"].add(read.query_name)

    bamfile_original_CT.close()
    

    bamfile_original_GA = pysam.AlignmentFile(bam_original_GA, 'rb')

    for read in bamfile_original_GA.fetch():
        hs = int(hashlib.sha224(read.query_name.encode()).hexdigest(), 16) % 6091

        if read.query_name not in dict_reads["GA"][hs]: continue

        if read.flag & 256 != 256 and read.flag & 2048 != 2048: continue

        if read.flag & 256 == 256:
            multi_reads["GA"].add(read.query_name)

        if read.flag & 2048 == 2048:
            split_reads["GA"].add(read.query_name)

    bamfile_original_GA.close()

    print('# of reads with multiple hits: ', len(multi_reads["CT"]) + len(multi_reads["GA"]) ,sep='\t', file=sys.stderr)
    print('# of reads with split hits: ', len(split_reads["CT"]) + len(split_reads["GA"]) ,sep='\t', file=sys.stderr)

    if rname:
        for name in multi_reads["CT"]:
            print(name, 'Multiple hits', 'CT',  sep='\t')
        for name in multi_reads["GA"]:
            print(name, 'Multiple hits', 'GA',  sep='\t')

        for name in split_reads["CT"]:
            print(name, 'Split hits', 'CT',  sep='\t')
        for name in split_reads["GA"]:
            print(name, 'Split hits', 'GA',  sep='\t')

    detect_unique_alignment_nanoEM(multi_reads, bam_CT, bam_GA) 


def detect_unique_alignment_wgs(multi_reads:dict, bam:str) -> None:

    bamfile = pysam.AlignmentFile(bam, 'rb')

    uniq_reads = set()
    uniq_reads_over20 = set()

    for read in bamfile.fetch():

        if read.flag & 4 == 4 or read.flag & 256 == 256: continue

        if read.query_name in multi_reads: continue

        uniq_reads.add(read.query_name)

        if read.mapping_quality >= 20:
            uniq_reads_over20.add(read.query_name)

    bamfile.close()

    print('# of reads with unique alignment: ', len(uniq_reads), sep='\t', file=sys.stderr)
    print('# of reads with unique alignment and MAPQ >= 20: ', len(uniq_reads_over20), sep='\t', file=sys.stderr)

def detect_unique_alignment_nanoEM(multi_reads:dict, bam_CT:str, bam_GA:str) -> None:

    print('Counting C-to-T converted reads...', file=sys.stderr)
    detect_unique_alignment_wgs(multi_reads['CT'], bam_CT)

    print('Counting G-to-A converted reads...', file=sys.stderr)
    detect_unique_alignment_wgs(multi_reads['GA'], bam_GA)


def main():
    parser = argparse.ArgumentParser(description='Count reads with multiple hits and split hits')
    subparsers = parser.add_subparsers()

    parser_wgs = subparsers.add_parser('wgs_count', help='see `wgs_count -h`')
    parser_wgs.add_argument('--bam', '-b', required=True, help='BAM file')
    parser_wgs.add_argument('--name', '-n', action='store_true', help='Write a read name in stdout') 
    parser_wgs.set_defaults(fn=detect_from_wgs)


    parser_nanoEM = subparsers.add_parser('nanoEM_count', help='see `nanoEM_count -h`')
    parser_nanoEM.add_argument('--CT', '-C', required=True, help='BAM file of CT converted reads')
    parser_nanoEM.add_argument('--GA', '-G', required=True, help='BAM file of GA converted reads') 
    parser_nanoEM.add_argument('--original_CT', '-T', required=True, help='BAM file of original CT converted reads')
    parser_nanoEM.add_argument('--original_GA', '-A', required=True, help='BAM file of original GA converted reads') 
    parser_nanoEM.add_argument('--name', '-n', action='store_true', help='Write a read name in stdout') 
    parser_nanoEM.set_defaults(fn=detect_from_nanoEM)

    args = parser.parse_args()

    if args.fn == detect_from_wgs:
        detect_from_wgs(args.bam, args.name)
    elif args.fn == detect_from_nanoEM:
        detect_from_nanoEM(args.CT, args.GA, args.original_CT, args.original_GA, args.name)

if __name__ == '__main__':
    main()


