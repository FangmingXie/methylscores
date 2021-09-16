#!/usr/bin/env python3
"""
Scores include:
mCH CH mCH_rate mCH_fmr mCH_fumr mCH_mhl mCH_umhl 
mCG CG mCG_rate mCG_fmr mCG_fumr mCG_mhl mCG_umhl

fmr: fully methylated rate (read level)
fumr: fully unmethylated rate (read level)
mhl: methylation haplotype load (Gou et al 2017)
umhl: unmethylation haplotype load (symmetric to mhl)
mcconc: DNA methylation concurance (Shi et al 2021)
"""

import numpy as np
from scipy import linalg
import re
import time


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

def calc_mhl_fast(strings, letter, verbose=False):
    """MHL and uMHL -- same as calc_mhl; improved speed > 5x
    """
    if not strings:
        return [np.nan, np.nan]

    tfi = time.time()
    
    letter_lower = letter.lower()
    letter_upper = letter.upper()
    if verbose:
        print("f1", time.time()-tfi)

    length_counts = np.bincount([len(string) for string in strings])
    dim = len(length_counts) 
    if verbose:
        print("f2", time.time()-tfi)

    # methylated
    substrings_meth = re.split(r'{}+'.format(letter_lower), letter_lower.join(strings)) # get all substring including empty ones
    length_counts_meth = np.bincount([len(string) for string in substrings_meth], minlength=dim)
    
    # unmethylated
    substrings_unmeth = re.split(r'{}+'.format(letter_upper), letter_upper.join(strings)) # get all substring including empty ones
    length_counts_unmeth = np.bincount([len(string) for string in substrings_unmeth], minlength=dim)
    if verbose:
        print("f3", time.time()-tfi)
    
    length_counts_all = np.vstack([
        length_counts, 
        length_counts_meth, 
        length_counts_unmeth, 
        ])[:, 1:]
    if verbose:
        print("f4", time.time()-tfi)
    
    # transform matrix
    dim = dim - 1 # first entry - 0 length removed
    trans_mat = np.flip(linalg.hankel((np.arange(dim)+1)[::-1]), axis=0) # lower triangular matrix
    weights = (np.arange(dim)+1)/(dim*(dim+1)/2)
    
    # apply transformation 
    length_counts_all = np.dot(length_counts_all,trans_mat)
    fracs = length_counts_all[1:]/length_counts_all[0]
    mhls = np.dot(fracs, weights)
    
    if verbose:
        print("f6", time.time()-tfi)

    return mhls

def calc_mcconc(strings, letter):
    """Shi et al. 2021
    mcconc: number of unmethylated sites present in partially methylated reads
    """
    letter_lower = letter.lower()
    letter_upper = letter.upper()
    
    conc_sites = np.sum([s.count(letter_lower) 
                         for s in strings if letter_upper in s])
            
    return conc_sites