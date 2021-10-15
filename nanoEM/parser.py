#! /usr/bin/env python

import argparse
from .run import *

def create_parser():

    parser = argparse.ArgumentParser(prog = "nanoEM")

    subparsers = parser.add_subparsers()

    ##########
    # convert_ref
    convert_ref = subparsers.add_parser("convert_ref",
                                help = "Convert C/Gs in reference genome to T/As")
    
    convert_ref.add_argument("--ref", '-r', type=str,
                                help = "Path to reference genome")

    convert_ref.set_defaults(func = convert_ref_main)
    ##########


    ##########
    # alignment
    alignment = subparsers.add_parser("alignment",
                                help = "Align converted reads to converted reference genome and choose best alignments")

    alignment.add_argument("--thread", "-t", type=int, default=4, help = "Thread number")
    alignment.add_argument("--ref", "-r", type=str, 
                                help = "Path to converted reference genome")
    alignment.add_argument("--fastq", "-f", type=str, 
                                help = "Path to fastq file")  

    alignment.set_defaults(func = best_align_main)
    ##########

    ##########
    # call-methylation
    call_methylation = subparsers.add_parser("call-methylation",
                                help = "Convert C/Gs in a fastq file to T/As")

    call_methylation.add_argument("--reference", "-ref", 'action=store_true', 
                                help = "Reference for runnning sambamba")
    call_methylation.add_argument("--thread", "-t", type=int, default=4, help = "Thread number")
    call_methylation.add_argument("--ref_fasta", "-f", type=str, 
                                help = "Reference fasta file")
    call_methylation.add_argument("--cpgs", "-c", type=str, 
                                help = "CpG sites in reference genome [bed format]")
    call_methylation.add_argument("--out_prefix", "-o", type=str, 
                                help = "Prefix of output directory")
    call_methylation.set_defaults(func = call_mathylation_main)
    ##########

    return parser