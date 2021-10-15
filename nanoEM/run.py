#! /usr/bin/env python

import os, subprocess, shutil
from .convert_ref import *
from .convert_reads import *
from .best_align import *
from .call_methylation import *

def convert_ref_main(args):

    if args.ref[-3:] == '.fa':
        hout = args.ref[:-3] + '_converted.fa.gz'
    elif args.ref[-6:] == '.fa.gz':
        hout = hout = args.ref[:-6] + '_converted.fa.gz'    
    conversion_ref(args.ref, hout)


def best_align_main(args):
    
    conversion_reads(args.fastq)

    prefix = args.fastq.split('.fa')[0][:-3]

    if shutil.which("minimap2") is None:
        print("Minimap2 does not exist", file=sys.stderr)
        sys.exit(1)
    
    if shutil.which("samtools") is None:
        print("Samtools does not exist", file=sys.stderr)
        sys.exit(1)
    
    command = 'minimap2 -t ' + args.thread + ' --split-prefix temp_sam1 -ax map-ont ' + args.ref + '_converted.fa.gz ' + \
                             prefix + '_CT.fa  --eqx | samtools view -b | samtools sort -@ ' + args.thread + ' -o 1.sorted.bam'
    subprocess.check_call(command, shell=True)

    command = 'samtools index 1.sorted.bam'
    subprocess.check_call(command, shell=True)

    command = 'minimap2 -t ' + args.thread + ' --split-prefix temp_sam1 -ax map-ont ' + args.ref + '_converted.fa.gz ' + \
                             prefix + '_GA.fa --eqx | samtools view -b | samtools sort -@ ' + args.thread + ' -o 2.sorted.bam'
    subprocess.check_call(command, shell=True)

    command = 'samtools index 2.sorted.bam'
    subprocess.check_call(command, shell=True)

    out = dict()

    extract_information("1.sorted.bam", 1, out)
    extract_information("2.sorted.bam", 2, out)

    compare(out, "1.sorted.bam", args.fastq)


def call_mathylation_main(args):
    if shutil.which("sambamba") is None:
        print("Sambamba does not exist", file=sys.stderr)
        sys.exit(1)
    

    if args.ref:
        command = 'sambamba mpileup output_CT.sorted.bam -L ' + args.cpgs + ' -o ' +  args.out_prefix + 'pileup_CT.tsv -t ' + \
                   args.thread + ' --samtools -f ' + args.ref_fasta
        subprocess.check_call(command, shell=True)

        command = 'sambamba mpileup output_GA.sorted.bam -L ' + args.cpgs + ' -o ' +  args.out_prefix + 'pileup_GA.tsv -t ' + \
                   args.thread + ' --samtools -f ' + args.ref_fasta
        subprocess.check_call(command, shell=True)
    
    else:
        command = 'sambamba mpileup output_CT.sorted.bam -L ' + args.cpgs + ' -o ' +  args.out_prefix + 'pileup_CT.tsv -t ' + args.thread
        subprocess.check_call(command, shell=True)

        command = 'sambamba mpileup output_GA.sorted.bam -L ' + args.cpgs + ' -o ' +  args.out_prefix + 'pileup_GA.tsv -t ' + args.thread
        subprocess.check_call(command, shell=True)
    
    fw = args.out_prefix + 'pileup_CT.tsv'
    re = args.out_prefix + 'pileup_GA.tsv'
    
    out = dict()
    
    call_methylation(fw, out, '+', args.ref, args.cpgs)
    call_methylation(re, out, '-', args.ref, args.cpgs)

    w = open(args.out_prefix + 'methylation_frequency.tsv', 'w')
    print('chromosome', 'position', 'methyl', 'unmethyl', 'freqeucy', sep='\t', file=w) 
    for key in out:
        c, p = key.split(':')
        freq = out[key]['methyl'] / (out[key]['methyl'] + out[key]['unmethyl']) 
        print(c, p, out[key]['methyl'], out[key]['unmethyl'], freq, sep='\t', file=w) 