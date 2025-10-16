import json
import yaml
import pandas as pd
import pytest
from pytest import approx

with open('groups.airr.yaml') as f:
    groups = yaml.safe_load(f)

groups = groups['RepertoireGroup']
call_types = ['c_call', 'd_call', 'j_call', 'v_call']
levels = ['allele', 'gene', 'subgroup']
modes = ['unique', 'exists', 'proportion']
productivity = [True, False]
processing_stage = 'igblast.makedb'

def test_sequence_frequency_sum_groups():
    '''
    Tests to ensure that the sequence_frequency variable sums to 1 for each group, each call_type, and each level for unique values.
    The assertion checks if the sum is approximately 1 (using pytest.approx()), or exactly 0. The 0 is included to ensure that
    any empty levels of the input 
    '''
    for group in groups:
        group_id = group['repertoire_group_id']
        for call_type in call_types:
            for level in levels:
                orig_df = pd.read_csv(f'{group_id}.{processing_stage}.group.{call_type}.tsv', delimiter='\t')
                for mode in modes:
                    for productive in productivity:
                        df = orig_df[(orig_df['mode']==mode) & (orig_df['productive']==productive) & (orig_df['level']==level)]
                        # print(f'{group_id}\t{call_type}\t{level}\t{productive}\t{df.loc[0,'sequence_frequency']}')
                        assert df['sequence_frequency'].sum() == approx(1) or df['sequence_frequency'].sum() == 0, f"Expected the sum of"+\
                            f"sequence_frequency to be equal to approx 1 or exactly 0. Instead: \n"+\
                            f"∑(sequence_frequency) = {df['sequence_frequency'].sum()} \n"+\
                            f"Variables: \n\tgroup_id:{group_id} \n\tcall_type:{call_type} \n\tlevel:{level} \n\tmode:{mode} \n\tproductive:{productive}"

def test_duplicate_frequency_sum_groups():
    '''
    Tests to ensure that the sequence_frequency variable sums to 1 for each group, each call_type, and each level for unique values.
    The assertion checks if the sum is approximately 1 (using pytest.approx()), or exactly 0. The 0 is included to ensure that
    any empty levels of the input 
    '''
    for group in groups:
        group_id = group['repertoire_group_id']
        for call_type in call_types:
            for level in levels:
                orig_df = pd.read_csv(f'{group_id}.{processing_stage}.group.{call_type}.tsv', delimiter='\t')
                for mode in modes:
                    for productive in productivity:
                        df = orig_df[(orig_df['mode']==mode) & (orig_df['productive']==productive) & (orig_df['level']==level)]
                        # print(f'{group_id}\t{call_type}\t{level}\t{productive}\t{df.loc[0,'duplicate_frequency']}')
                        assert df['duplicate_frequency'].sum() == approx(1) or df['duplicate_frequency'].sum() == 0, f"Expected the sum of"+\
                            f"duplicate_frequency to be equal to approx 1 or exactly 0. Instead: \n"+\
                            f"∑(duplicate_frequency) = {df['duplicate_frequency'].sum()} \n"+\
                            f"Variables: \n\tgroup_id:{group_id} \n\tcall_type:{call_type} \n\tlevel:{level} \n\tmode:{mode} \n\tproductive:{productive}"
