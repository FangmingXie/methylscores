#!/usr/bin/env python3
"""
Scores include:
mCH CH mCH_rate mCH_fmr mCH_fumr mCH_mhl mCH_umhl 
mCG CG mCG_rate mCG_fmr mCG_fumr mCG_mhl mCG_umhl

fmr: fully methylated rate (read level)
fumr: fully unmethylated rate (read level)
mhl: methylation haplotype load (Gou et al 2017)
umhl: unmethylation haplotype load (symmetric to mhl)
"""

import numpy as np
from scipy import linalg
import pandas as pd
import argparse
import re
import logging
import pysam
import time
from multiprocessing import Pool

import compute_scores_utils # current directory


# prototype6 -  wrap prototype5 up and allow multiprocessing
def calculate_read_level_mc_metrics_genomewide(
    input_bam, output_scores, 
    bin_size=1000, # 1kb
    verbose_level=1000000 # update every 1 million reads
    ):
    
    """
    Requirements:
    - sorted bam file
    - bam file created by Bismark (has the "XM" tag in the third position of tags)
    - bin_size > sequence length
    """
    logging.info("input: {}".format(input_bam))
    logging.info("output: {}".format(output_scores))
    
    ti = time.time()
    columns = [
        'chr', 'start', 'end',
        'ch_mc', 'ch_c',
        'cg_mc', 'cg_c',
        'ch_fully_meth_reads', 'ch_fully_unmeth_reads', 'ch_total_reads',
        'cg_fully_meth_reads', 'cg_fully_unmeth_reads', 'cg_total_reads',
        'ch_mhl', 'ch_umhl',
        'cg_mhl', 'cg_umhl',
    ]

    bamfile = pysam.AlignmentFile(input_bam, 'rb')
    assert bamfile.has_index()

    current_chrom = "chr0"
    current_binstart = 0
    current_binend = current_binstart + bin_size 
    current_binstring = ""
    next_binstring_hold = "" 

    with open(output_scores, 'w') as output_fh:
        # title
        output_fh.write("\t".join(columns)+'\n')
        for i, read in enumerate(bamfile.fetch()):
            if i % verbose_level == 0:
                logging.info("{} {} {}".format(input_bam.split('/')[-1], i, time.time()-ti))

            # deal with one read
            chrom, start, end, mc_string = (read.reference_name, 
                                            read.reference_start, 
                                            read.reference_end, 
                                            read.tags[2][1],
                                           )
            # handle chr
            if not chrom.startswith('chr'):
                chrom = 'chr'+chrom

            # when current read goes beyond the current bin
            # keep updating the current bin until the current read is within the bin 
            while (current_chrom != chrom) or (start > current_binend):
                # flush the current bin out
                if current_binstring:
                    # compute all scores
                    mch_str, mcg_str = compute_scores_utils.process_strings(current_binstring)
                    # site level
                    ch_mc, ch_c = compute_scores_utils.calc_mc(mch_str, 'h')
                    cg_mc, cg_c = compute_scores_utils.calc_mc(mcg_str, 'z')
                    # read level
                    ch_fmr, ch_fumr, ch_total_reads = compute_scores_utils.calc_fmr(mch_str, 'h')
                    cg_fmr, cg_fumr, cg_total_reads = compute_scores_utils.calc_fmr(mcg_str, 'z')
                    # MHL level
                    ch_mhl, ch_umhl = compute_scores_utils.calc_mhl_fast(mch_str, 'h')
                    cg_mhl, cg_umhl = compute_scores_utils.calc_mhl_fast(mcg_str, 'z')

                    # record values
                    column_values = [
                        current_chrom, current_binstart, current_binend,
                        ch_mc, ch_c,
                        cg_mc, cg_c,
                        ch_fmr, ch_fumr, ch_total_reads,
                        cg_fmr, cg_fumr, cg_total_reads,
                        ch_mhl, ch_umhl,
                        cg_mhl, cg_umhl,
                    ]
                    # output 
                    output_fh.write("\t".join([str(item) for item in column_values])+'\n')

                # move current bin to the next 
                # (move to a new chrom or the next bin within the chrom)
                if current_chrom != chrom:
                    current_chrom = chrom
                    current_binstart = 0 
                    current_binend = bin_size
                    current_binstring = next_binstring_hold
                    next_binstring_hold = ""
                else:
                    current_binstart += bin_size
                    current_binend += bin_size
                    current_binstring = next_binstring_hold
                    next_binstring_hold = ""

            # this means current_chrom == chrom, and start <= current bin
            if end <= current_binend:
                # move the whole read into the current bin
                current_binstring += (mc_string+',')
            else:
                # move the first section of the read into the current bin
                current_binstring += (mc_string[:-(end-current_binend)]+',')
                # hold the second section (beyond the current bin) temporarily
                next_binstring_hold += (mc_string[(current_binend-start):]+',')

        # looped over all reads
    bamfile.close()
    logging.info("Total time: {:.2f} ({})".format(time.time()-ti, input_bam))
    
    return

def create_parser():
    """
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--inputs", 
    	required=True,
    	nargs="+",
    	help="sorted bam file generated from Bismark")
    parser.add_argument("-o", "--outputs", 
    	required=True,
    	nargs="+",
    	help="bed format")
    parser.add_argument("-s", "--bin_size", 
    	type=int,
    	required=True,
    	help="bin size")
    parser.add_argument("-v", "--verbose_level", 
    	type=int,
    	default=1000000,
    	help="update every (1 million) reads")
    parser.add_argument("-n", "--ncores", 
    	type=int,
    	default=1,
    	help="number of cores to use")

    return parser

def create_logger(name='log'):
    """
    args: logger name

    return: a logger object
    """
    logging.basicConfig(
        format='%(asctime)s %(message)s', 
        datefmt='%m/%d/%Y %I:%M:%S %p',
        level=logging.INFO)
    return logging.getLogger(name)


if __name__ == "__main__":

	log = create_logger()
	parser = create_parser()
	args = parser.parse_args()

	inputs = args.inputs
	outputs = args.outputs
	bin_size = args.bin_size
	verbose_level = args.verbose_level
	ncores = args.ncores

	logger = create_logger()
	logging.info("Calculate mc metrics genomewide")
	logging.info("{} inputs, bin size {}, verbose_level {}, ncores {}".format(
				len(inputs), bin_size, verbose_level, ncores))


	with Pool(min(ncores, len(inputs))) as p:
	#     # test 1 
	#     calculate_read_level_mc_metrics_genomewide(
	#                         _input, _output, 
	#                         bin_size=bin_size,
	#                         verbose_level=verbose_level,
	#                         )
	    # parallel
	    p.starmap(calculate_read_level_mc_metrics_genomewide, 
	                  [(_input, _output, bin_size, verbose_level,) 
	                       for _input, _output in zip(inputs, outputs)],
	             )