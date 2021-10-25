#! /usr/bin/env python

import os, subprocess, shutil
from convert_ref import *
from convert_reads import *
from best_align import *
from call_methylation import *

def convert_ref_main(args):
    
    conversion_ref(args.ref, args.out_prefix, args.thread)


def best_align_main(args):
    
    conversion_reads(args.fastq)

    prefix = args.fastq.split('.fq')[0]

    if shutil.which("minimap2") is None:
        print("Minimap2 does not exist", file=sys.stderr)
        sys.exit(1)
    
    if shutil.which("samtools") is None:
        print("Samtools does not exist", file=sys.stderr)
        sys.exit(1)
    
    if not os.path.exists(args.out_prefix):
        print("Output directory does not exit!", file=sys.stderr)
        sys.exit(1)

    command = 'minimap2 -t ' + str(args.thread) + ' --split-prefix ' + args.out_prefix + 'temp_sam1 -ax map-ont ' + args.ref + ' ' + \
                             prefix + '_CT.fq  --eqx | samtools view -b | samtools sort -@ ' + str(args.thread) + ' -o ' + args.out_prefix + '1.sorted.bam'
    subprocess.check_call(command, shell=True)

    command = 'samtools index ' + args.out_prefix + '1.sorted.bam'
    subprocess.check_call(command, shell=True)

    command = 'minimap2 -t ' + str(args.thread) + ' --split-prefix ' + args.out_prefix + 'temp_sam1 -ax map-ont ' + args.ref + ' ' + \
                             prefix + '_GA.fq --eqx | samtools view -b | samtools sort -@ ' + str(args.thread) + ' -o ' + args.out_prefix + '2.sorted.bam'
    subprocess.check_call(command, shell=True)

    command = 'samtools index ' + args.out_prefix + '2.sorted.bam'
    subprocess.check_call(command, shell=True)

    out = dict()

    extract_information(args.out_prefix + "1.sorted.bam", 1, out)
    extract_information(args.out_prefix + "2.sorted.bam", 2, out)

    compare(out, args.out_prefix + "1.sorted.bam", args.fastq, args.out_prefix, args.thread)

    os.remove(args.out_prefix + '1.sorted.bam')
    os.remove(args.out_prefix + '1.sorted.bam.bai')
    os.remove(args.out_prefix + '2.sorted.bam')
    os.remove(args.out_prefix + '2.sorted.bam.bai')


def call_mathylation_main(args):
    if shutil.which("sambamba") is None:
        print("Sambamba does not exist", file=sys.stderr)
        sys.exit(1)
    

    if args.reference:
        #command = 'sambamba mpileup output_CT.sorted.bam -L ' + args.cpgs + ' -o ' +  args.out_prefix + 'pileup_CT.tsv -t ' + \
        #           args.thread + ' --samtools -f ' + args.ref_fasta
        #subprocess.check_call(command, shell=True)

        hout = open(args.out_prefix + 'pileup_CT.tsv', 'w')
        subprocess.check_call(['samtools', 'mpileup', '-l', args.cpgs, '-f', args.ref_fasta, args.input_prefix + 'output_CT.sorted.bam'], stdout = hout)
        hout.close()

        #command = 'sambamba mpileup output_GA.sorted.bam -L ' + args.cpgs + ' -o ' +  args.out_prefix + 'pileup_GA.tsv -t ' + \
        #           args.thread + ' --samtools -f ' + args.ref_fasta
        #subprocess.check_call(command, shell=True)
        hout = open(args.out_prefix + 'pileup_GA.tsv', 'w')
        subprocess.check_call(['samtools', 'mpileup', '-l', args.cpgs, '-f', args.ref_fasta, args.input_prefix + 'output_GA.sorted.bam'], stdout = hout)
        hout.close()
    
    else:
        #command = 'sambamba mpileup output_CT.sorted.bam -L ' + args.cpgs + ' -o ' +  args.out_prefix + 'pileup_CT.tsv -t ' + args.thread
        #subprocess.check_call(command, shell=True)
        hout = open(args.out_prefix + 'pileup_CT.tsv', 'w')
        subprocess.check_call(['samtools', 'mpileup', '-l', args.cpgs, 'output_CT.sorted.bam'], stdout = hout)
        hout.close()

        #command = 'sambamba mpileup output_GA.sorted.bam -L ' + args.cpgs + ' -o ' +  args.out_prefix + 'pileup_GA.tsv -t ' + args.thread
        #subprocess.check_call(command, shell=True)

        hout = open(args.out_prefix + 'pileup_GA.tsv', 'w')
        subprocess.check_call(['samtools', 'mpileup', '-l', args.cpgs, 'output_GA.sorted.bam'], stdout = hout)
        hout.close()
    
    fw = args.out_prefix + 'pileup_CT.tsv'
    re = args.out_prefix + 'pileup_GA.tsv'
    
    out = dict()
    
    call_methylation(fw, out, '+', args.reference, args.cpgs)
    call_methylation(re, out, '-', args.reference, args.cpgs)

    w = open(args.out_prefix + 'methylation_frequency.tsv', 'w')
    print('chromosome', 'position', 'methyl', 'unmethyl', 'freqeucy', sep='\t', file=w) 
    for key in out:
        c, p = key.split(':')
        if out[key]['methyl'] + out[key]['unmethyl']: continue
        freq = out[key]['methyl'] / (out[key]['methyl'] + out[key]['unmethyl']) 
        print(c, p, out[key]['methyl'], out[key]['unmethyl'], freq, sep='\t', file=w) 