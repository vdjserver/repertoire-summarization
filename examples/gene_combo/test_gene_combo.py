import os
import json
import yaml
import pandas as pd
import pytest
from pytest import approx

curr_dir = os.path.dirname(__file__)
with open(os.path.join(curr_dir, 'groups.airr.yaml')) as f:
    groups = yaml.safe_load(f)

groups = groups['RepertoireGroup']
combos = ['dj_combo', 'vd_combo', 'vdj_combo', 'vj_combo']
levels = ['allele|allele', 'gene|gene', 'subgroup|subgroup']
modes = ['unique', 'exists', 'proportion']
productivity = [True, False]
processing_stage = 'igblast.makedb'

def test_sequence_frequency_sum_groups():
    '''
    Tests to ensure that the sequence_frequency variable sums to 1 for each group, each call_type, and each level for unique values.
    The assertion checks if the sum is approximately 1 (using pytest.approx()), or exactly 0.
    '''
    for group in groups:
        group_id = group['repertoire_group_id']
        rep_ids = [rep_obj['repertoire_id'] for rep_obj in group['repertoires']]
        for combo in combos:
            for level in levels:
                filename = f'{group_id}.{processing_stage}.group.repertoires.{combo}.tsv'
                orig_df = pd.read_csv(os.path.join(curr_dir, filename), delimiter='\t')
                for mode in modes:
                    for productive in productivity:
                        for rep_id in rep_ids:
                            df = orig_df[(orig_df['mode']==mode) & (orig_df['productive']==productive) & (orig_df['level']==level) & (orig_df['repertoire_id']==rep_id)]
                        
                            assert df['sequence_frequency'].sum() == approx(1) or df.empty, f"Expected the sum of"+\
                                f"sequence_frequency to be equal to approx 1 or for dataframe to be empty. Instead: \n"+\
                                f"∑(sequence_frequency) = {df['sequence_frequency'].sum()} \n"+\
                                f"Variables: \n\tgroup_id:{group_id} \n\trepertoire_id:{rep_id} \n\tcall_type:{combo} \n\tlevel:{level} \n\tmode:{mode} \n\tproductive:{productive}"

def test_duplicate_frequency_sum_groups():
    '''
    Tests to ensure that the sequence_frequency variable sums to 1 for each group, each call_type, and each level for unique values.
    The assertion checks if the sum is approximately 1 (using pytest.approx()), or for the dataframe to be empty.
    '''
    for group in groups:
        group_id = group['repertoire_group_id']
        rep_ids = [rep_obj['repertoire_id'] for rep_obj in group['repertoires']]
        for combo in combos:
            for level in levels:
                filename = f'{group_id}.{processing_stage}.group.repertoires.{combo}.tsv'
                orig_df = pd.read_csv(os.path.join(curr_dir, filename), delimiter='\t')
                for mode in modes:
                    for productive in productivity:
                        for rep_id in rep_ids:
                            df = orig_df[(orig_df['mode']==mode) & (orig_df['productive']==productive) & (orig_df['level']==level) & (orig_df['repertoire_id']==rep_id)]
                        
                            assert df['duplicate_frequency'].sum() == approx(1) or df.empty, f"Expected the sum of"+\
                                f"duplicate_frequency to be equal to approx 1 or for dataframe to be empty. Instead: \n"+\
                                f"∑(duplicate_frequency) = {df['duplicate_frequency'].sum()} \n"+\
                                f"Variables: \n\tgroup_id:{group_id} \n\trepertoire_id:{rep_id} \n\tcall_type:{combo} \n\tlevel:{level} \n\tmode:{mode} \n\tproductive:{productive}"
