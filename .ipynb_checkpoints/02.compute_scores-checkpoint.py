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
import logging
import re


def process_strings(strings, consider_pairedend=True):
    """Process strings
    split mCH and mCG
    z Z - mCG
    x X - mCHG, h H - mCHH
    u U - mCN
    . not C
    , read deliminator
    ; paired end read deliminator

    consider_pairedend - True: combine paired end reads as 1 fragment
                       - False: consider them as 2 separate single-end reads 
                       
    return: hH list/zZ list
    """
    if consider_pairedend:
        strings = (strings.replace('.', '')
                          .replace('u', '')
                          .replace('x', 'h')
                          .replace('X', 'H')
                          .replace(';', '') # consider paired-end reads are linked
                   )
    else:
        strings = (strings.replace('.', '')
                          .replace('u', '')
                          .replace('x', 'h')
                          .replace('X', 'H')
                          .replace(';', ',') # consider paired-end reads are separate 
                   )

    # z Z h H
    string_list_mch = (strings.replace('z', '')
                             .replace('Z', '')
                             .split(',')
                      )
    string_list_mcg = (strings.replace('h', '')
                             .replace('H', '')
                             .split(',')
                      )
    
    # remove empty entries
    string_list_mch = [string for string in string_list_mch if string]
    string_list_mcg = [string for string in string_list_mcg if string]
    return string_list_mch, string_list_mcg

def calc_mc(string_list, letter):
    """
    """
    letter = letter.upper()
    string = ''.join(string_list)
    mc = string.count(letter)
    c = len(string)
    return mc, c

def calc_fmr(string_list, letter):
    """
    """
    letter_upper = letter.upper()
    letter_lower = letter.lower()
    string_sets = np.array([set(string) for string in string_list])
    num_fm = (string_sets == {letter_upper}).sum()
    num_fum = (string_sets == {letter_lower}).sum()
    num = len(string_sets)
    
    return num_fm, num_fum, num 

def calc_mhl(strings, letter):
    """MHL and uMHL
    """
    letter_lower = letter.lower()
    letter_upper = letter.upper()

    length_counts = np.bincount([len(string) for string in strings])
    length_counts_meth = np.zeros_like(length_counts)
    length_counts_unmeth = np.zeros_like(length_counts)

    for string in strings:
        # methylated
        a = re.split(r'{}+'.format(letter_lower), string.strip(letter_lower)) # get all substring including empty ones
        a = np.bincount([len(_a) for _a in a])
        length_counts_meth[:len(a)] += a

        # unmethylated
        a = re.split(r'{}+'.format(letter_upper), string.strip(letter_upper)) # get all substring including empty ones
        a = np.bincount([len(_a) for _a in a])
        length_counts_unmeth[:len(a)] += a

    length_counts_all = np.vstack([
        length_counts, 
        length_counts_meth, 
        length_counts_unmeth, 
        ])[:, 1:]

    dim = length_counts_all.shape[1]
    trans_mat = np.flip(linalg.hankel((np.arange(dim)+1)[::-1]), axis=0) # lower triangular matrix
    length_counts_all = np.dot(length_counts_all,trans_mat)

    fracs = length_counts_all[1:]/length_counts_all[0]
    weights = (np.arange(dim)+1)/(dim*(dim+1)/2)
    mhls = np.dot(fracs, weights)

    return mhls

def main(input_mcinfo, output_mcscores):
	"""Calc all scores
	"""
	df = pd.read_csv(input_mcinfo, sep='\t', header=None,
	                 names=['chr', 'start', 'end', 'seq'], 
	                )
	logging.info("input size: {}".format(df.shape))

	df = df[~df['seq'].isnull()]
	logging.info("input size (net): {}".format(df.shape))



	logging.info("Process reads info...")
	df_proc = pd.DataFrame(df[~df['seq'].isnull()]['seq'].apply(process_strings).tolist(), 
	                       index=df['seq'].index,
	                       columns=['mch_string', 'mcg_string'])

	logging.info("Calculate mC at site level...")
	df_mch = pd.DataFrame(df_proc['mch_string'].apply(lambda x: calc_mc(x, 'h')).tolist(),
	                     index=df_proc['mch_string'].index,
	                     columns=['ch_mc', 'ch_c'],
	                     )
	df_mcg = pd.DataFrame(df_proc['mcg_string'].apply(lambda x: calc_mc(x, 'z')).tolist(),
	                     index=df_proc['mcg_string'].index,
	                     columns=['cg_mc', 'cg_c'],
	                     )

	logging.info("Calculate mC at read level...")
	df_mch_fmr = pd.DataFrame(df_proc['mch_string'].apply(lambda x: calc_fmr(x, 'h')).tolist(),
	                     index=df_proc['mch_string'].index,
	                     columns=['ch_fully_meth_reads', 'ch_fully_unmeth_reads', 'ch_total_reads'],
	                    )

	df_mcg_fmr = pd.DataFrame(df_proc['mcg_string'].apply(lambda x: calc_fmr(x, 'z')).tolist(),
	                     index=df_proc['mcg_string'].index,
	                     columns=['cg_fully_meth_reads', 'cg_fully_unmeth_reads', 'cg_total_reads'],
	                    )

	logging.info("Calculate MHL...")
	df_mch_mhl = pd.DataFrame(df_proc['mch_string'].apply(lambda x: calc_mhl(x, 'h')).tolist(),
	                     index=df_proc['mch_string'].index,
	                     columns=['ch_mhl', 'ch_umhl'],
	                    )
	df_mcg_mhl = pd.DataFrame(df_proc['mcg_string'].apply(lambda x: calc_mhl(x, 'z')).tolist(),
	                     index=df_proc['mcg_string'].index,
	                     columns=['cg_mhl', 'cg_umhl'],
	                    )

	logging.info("Gather results...")
	df_res = pd.concat([ 
	                   df.drop('seq', axis=1),
	                   df_mch, df_mcg, 
	                   df_mch_fmr, df_mcg_fmr, 
	                   df_mch_mhl, df_mcg_mhl,
	                  ], axis=1)

	logging.info("Save results...")
	df_res.to_csv(output_mcscores, 
					sep='\t', 
					header=True, index=True, 
					na_rep='NA')
 
def create_parser():
    """
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_mcinfo", 
    	required=True,
    	help="bed format")
    parser.add_argument("-o", "--output_mcscores", 
    	required=True,
    	help="bed format")
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

if __name__ == '__main__':

	log = create_logger()
	parser = create_parser()
	args = parser.parse_args()

	input_mcinfo = args.input_mcinfo
	output_mcscores = args.output_mcscores
	logging.info("input: {}\noutput: {}"
					.format(input_mcinfo, output_mcscores))

	main(input_mcinfo, output_mcscores)
