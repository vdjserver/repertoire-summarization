'''
Author: Sam Wollenburg
Date: Mar. 18, 2025

To be used after repcalc mutations.py has been run. The input of this file is the output from mutations.py
This calculates various metrics regarding RGWY/WRCY (RW) motifs. Metric captured include: Total number of 
RW motifs in sequence, total number of RW motifs contained within a FWR or CDR region, the number of mutated
nucleotides in the sequence and the number of mutated nucleotides for RW motifs contained within a FWR or CDR
region. Calculates frequency of mutated nucleotides in RW motifs in sequence and FWR and CDR regions.
'''

from repcalc import __version__
import repcalc.defaults as defaults
import repcalc.metadata as metadata
# import repcalc.gldb as gldb
import json
import math
import numpy as np
import csv
import airr
from Bio.Seq import Seq
import pandas as pd

# calculate post-`mutations.py`
samples = pd.read_csv('/Users/s236922/code/projects/analysis_v3/gt3shm/UTSW33_S42_L001_R1_001.fastq.merged.unique.gt3shm.mutations.airr.tsv', delimiter='\t')
# reduce size for necessary components to use and for write out
samples_count = samples.iloc[:,0:9]
samples_freq = samples.iloc[:,0:9]

# [start nt idx, end nt idx) 
regions = {'fwr1':(0, 78), 'cdr1':(78, 117), 'fwr2':(117, 168), 'cdr2':(168, 198), 'fwr3':(198, 315)}

# rw motifs
rgwy = ['AGCT', 'AGCA', 'AGTT', 'AGTA', 'GGCT', 'GGCA', 'GGTT', 'GGTA']
wrcy = ['TACC', 'TACT', 'TGCC', 'TGCT', 'AACC', 'AACT', 'AGCC', 'AGCT']

# new count columns to add
new_count_cols = ['mu_total_count_rw_motif', 'mu_count_rw_motif',
                  'mu_total_count_rw_motif_fwr1', 'mu_count_rw_motif_fwr1', 'mu_total_count_rw_motif_cdr1', 'mu_count_rw_motif_cdr1',
                  'mu_total_count_rw_motif_fwr2', 'mu_count_rw_motif_fwr2', 'mu_total_count_rw_motif_cdr2', 'mu_count_rw_motif_cdr2',
                  'mu_total_count_rw_motif_fwr3', 'mu_count_rw_motif_fwr3']

# new frequency columns to add
new_freq_cols = ['mu_freq_rw_motif',
                 'mu_freq_rw_motif_fwr1', 'mu_freq_rw_motif_cdr1',
                 'mu_freq_rw_motif_fwr2', 'mu_freq_rw_motif_cdr2',
                 'mu_freq_rw_motif_fwr3']

# init new columns with `None`
samples_count[new_count_cols] = np.full([samples.shape[0], len(new_count_cols)], None)
samples_freq[new_freq_cols] = np.full([samples.shape[0], len(new_freq_cols)], None) 

def find_rw_motifs(sequence:str, init_cutoff:int=45) -> dict:
    '''
    Finds the numer of rgwy and wrcy motifs in a sequence

    Parameters
    ----------
    sequence : str
        from germline
    init_cutoff : int
        the amount of unavailable sequence to remove from beginning

    Returns
    -------
    rw_idx : dict
        Dictionary with idx as key, motif as value of the rw motif index locations
        
    '''
     # amount of sequence that is removed due to not being available
    rw_idx = {}

    # loop through all 4-mer in sequence from init_cutoff
    for i in range(init_cutoff, len(sequence)-3):
        
        # Ensure motif doesn't start with a gap
        if sequence[i] != '.':
            # extract a sequence that is 4 nt long w/ gaps
            i_end = i + 3
            nt_count = 0
            while nt_count < 4 and i_end < len(sequence)-1:
                # remove gap for counting
                i_end += 1
                kmer = ''.join(sequence[i:i_end].split('.'))
                nt_count = len(kmer)
                

            # compare gapless w/ list of motifs
            if kmer in rgwy or kmer in wrcy:
                # add gap sequence to dictionary
                rw_idx[i] = sequence[i:i_end]
    
    return rw_idx

def find_rw_mutations(sequence:str, rw_idx:dict) -> dict:
    '''
    Finds the number of rw motif mutations from a sample's sequence, and a list of know rw motif indicies from germline.

    Parameters
    ----------
    sequence : str
        from sample
    rw_idx : dict
        from find_rw_motifs(). {motif_index : motif}
    
    Returns
    -------
    rw_mutations : dict
        format {'region_name':mutation_count} {string: int}
    '''
    rw_mutations = {'all':0,
                    'fwr1':0, 'cdr1':0,
                    'fwr2':0, 'cdr2':0,
                    'fwr3':0}
    
    def inner_func(curr_region='all'):
        # loop through all rw locations
        for i_start, motif_gl in rw_idx.items():
            # if different look at motif to see # of mutations
            if sequence[i_start:i_start+len(motif_gl)] != motif_gl:
                # across whole sequence
                if curr_region == 'all':
                    # march through sequence and count mutations
                    for i in range(len(motif_gl)):
                        if sequence[i_start:i_start+len(motif_gl)][i] != motif_gl[i]:
                            rw_mutations[curr_region] += 1
                # across regions
                else:
                    # if sequence is within a region 
                    if i_start >= regions[curr_region][0] and i_start < regions[curr_region][1]-len(motif_gl):
                        # march through sequence and count mutations
                        for i in range(len(motif_gl)):
                            if sequence[i_start:i_start+len(motif_gl)][i] != motif_gl[i]:
                                rw_mutations[curr_region] += 1
                        

    
    # find all rw motif mutation count in sequence
    inner_func('all')

    # find each regions rw motif mutation count
    for region_name in regions.keys():
        inner_func(region_name)
    
    return rw_mutations

def find_rw_distribution(rw_idx:dict)->dict:
    '''
    Tallies-up the number of rw motifs contained within a region 
    Ignores overlap between regions
 
    Parameters
    ----------
    rw_idx : dict
        from find_rw_motifs(). {motif_index : motif}
    
    Returns
    -------
    distribution : dict
        contains the distribution as {<region name>:<rw motif count>}
    '''
    distribution = {'fwr1':0, 'cdr1':0,
                    'fwr2':0, 'cdr2':0,
                    'fwr3':0}
    
    # Move through each region
    for region_name, region_start_end in regions.items():
        for idx, rw in rw_idx.items():
        
            if idx >= region_start_end[0] and idx < region_start_end[1]-len(rw):
                distribution[region_name] += 1
    
    return distribution

def main():
    # update count dataframe
    for sample in samples_count.itertuples():
        # create dict of rw motif idx in germline
        rw_idx = find_rw_motifs(sample.germline_alignment)

        # count motif breakdown by region
        motif_distrib = find_rw_distribution(rw_idx)

        # find mutation count in each region
        rw_mutations = find_rw_mutations(sample.sequence_alignment, rw_idx)
        
        # count number of rw motifs in germline
        samples_count.loc[sample.Index, 'mu_total_count_rw_motif'] = len(rw_idx)
        samples_count.loc[sample.Index, 'mu_total_count_rw_motif_fwr1'] = motif_distrib['fwr1']
        samples_count.loc[sample.Index, 'mu_total_count_rw_motif_cdr1'] = motif_distrib['cdr1']
        samples_count.loc[sample.Index, 'mu_total_count_rw_motif_fwr2'] = motif_distrib['fwr2']
        samples_count.loc[sample.Index, 'mu_total_count_rw_motif_cdr2'] = motif_distrib['cdr2']
        samples_count.loc[sample.Index, 'mu_total_count_rw_motif_fwr3'] = motif_distrib['fwr3']

        # count number of rw mutations in sample compared to germline
        samples_count.loc[sample.Index, 'mu_count_rw_motif'] = rw_mutations['all']
        samples_count.loc[sample.Index, 'mu_count_rw_motif_cdr1'] = rw_mutations['fwr1']
        samples_count.loc[sample.Index, 'mu_count_rw_motif_fwr1'] = rw_mutations['cdr1']
        samples_count.loc[sample.Index, 'mu_count_rw_motif_fwr2'] = rw_mutations['fwr2']
        samples_count.loc[sample.Index, 'mu_count_rw_motif_cdr2'] = rw_mutations['cdr2']
        samples_count.loc[sample.Index, 'mu_count_rw_motif_fwr3'] = rw_mutations['fwr3']

    # update frequence dataframe
    for sample in samples_freq.itertuples():
        if samples_count.loc[sample.Index, 'mu_total_count_rw_motif'] != 0:
            samples_freq.loc[sample.Index, 'mu_freq_rw_motif'] = samples_count.loc[sample.Index, 'mu_count_rw_motif']/(samples_count.loc[sample.Index, 'mu_total_count_rw_motif']*4)
        else:
            samples_freq.loc[sample.Index, 'mu_freq_rw_motif'] = np.NaN

        if samples_count.loc[sample.Index, 'mu_total_count_rw_motif_fwr1'] != 0:
            samples_freq.loc[sample.Index, 'mu_freq_rw_motif_fwr1'] = samples_count.loc[sample.Index, 'mu_count_rw_motif_fwr1']/(samples_count.loc[sample.Index, 'mu_total_count_rw_motif_fwr1']*4)
        else:
            samples_freq.loc[sample.Index, 'mu_freq_rw_motif_fwr1'] = np.NaN

        if samples_count.loc[sample.Index, 'mu_total_count_rw_motif_cdr1'] != 0:
            samples_freq.loc[sample.Index, 'mu_freq_rw_motif_cdr1'] = samples_count.loc[sample.Index, 'mu_count_rw_motif_cdr1']/(samples_count.loc[sample.Index, 'mu_total_count_rw_motif_cdr1']*4)
        else:
            samples_freq.loc[sample.Index, 'mu_freq_rw_motif_cdr1'] = np.NaN

        if samples_count.loc[sample.Index, 'mu_total_count_rw_motif_fwr2'] != 0:
            samples_freq.loc[sample.Index, 'mu_freq_rw_motif_fwr2'] = samples_count.loc[sample.Index, 'mu_count_rw_motif_fwr2']/(samples_count.loc[sample.Index, 'mu_total_count_rw_motif_fwr2']*4)
        else:
            samples_freq.loc[sample.Index, 'mu_freq_rw_motif_fwr2'] = np.NaN

        if samples_count.loc[sample.Index, 'mu_total_count_rw_motif_cdr2'] != 0:
            samples_freq.loc[sample.Index, 'mu_freq_rw_motif_cdr2'] = samples_count.loc[sample.Index, 'mu_count_rw_motif_cdr2']/(samples_count.loc[sample.Index, 'mu_total_count_rw_motif_cdr2']*4)
        else:
            samples_freq.loc[sample.Index, 'mu_freq_rw_motif_cdr2'] = np.NaN

        if samples_count.loc[sample.Index, 'mu_total_count_rw_motif_fwr3'] != 0:
            samples_freq.loc[sample.Index, 'mu_freq_rw_motif_fwr3'] = samples_count.loc[sample.Index, 'mu_count_rw_motif_fwr3']/(samples_count.loc[sample.Index, 'mu_total_count_rw_motif_fwr3']*4)
        else:
            samples_freq.loc[sample.Index, 'mu_freq_rw_motif_fwr3'] = np.NaN
    
    samples_count.to_csv('/Users/s236922/code/projects/analysis_v3/gt3shm/UTSW33_S42_L001_R1_001.fastq.merged.unique.gt3shm.mutations.rw_motif_count.airr.tsv', delimiter='\t')
    samples_freq.to_csv('/Users/s236922/code/projects/analysis_v3/gt3shm/UTSW33_S42_L001_R1_001.fastq.merged.unique.gt3shm.mutations.rw_motif_freq.airr.tsv', delimiter='\t')

