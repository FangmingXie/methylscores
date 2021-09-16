#!/usr/bin/env python3
"""Wrappers around compute_scores_utils
"""

import numpy as np
import pandas as pd
import argparse
import logging

from compute_scores_utils import *


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
    df_mch_mhl = pd.DataFrame(df_proc['mch_string'].apply(lambda x: calc_mhl_fast(x, 'h')).tolist(),
                         index=df_proc['mch_string'].index,
                         columns=['ch_mhl', 'ch_umhl'],
                        )
    df_mcg_mhl = pd.DataFrame(df_proc['mcg_string'].apply(lambda x: calc_mhl_fast(x, 'z')).tolist(),
                         index=df_proc['mcg_string'].index,
                         columns=['cg_mhl', 'cg_umhl'],
                        )

    logging.info("Calculate mC concurance...")
    df_mch_conc = pd.DataFrame(df_proc['mch_string'].apply(lambda x: calc_mcconc(x, 'h')).tolist(),
                         index=df_proc['mch_string'].index,
                         columns=['ch_conc'],
                        )
    df_mcg_conc = pd.DataFrame(df_proc['mcg_string'].apply(lambda x: calc_mcconc(x, 'z')).tolist(),
                         index=df_proc['mcg_string'].index,
                         columns=['cg_conc'],
                        )

    logging.info("Gather results...")
    df_res = pd.concat([ 
                       df.drop('seq', axis=1),
                       df_mch, df_mcg, 
                       df_mch_fmr, df_mcg_fmr, 
                       df_mch_mhl, df_mcg_mhl,
                       df_mch_conc, df_mcg_conc,
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
