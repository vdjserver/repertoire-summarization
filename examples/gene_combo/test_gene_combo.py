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

def test_field_sum_group_repertoires(field):
    '''
    Tests to ensure that the given field variable sums to 1 for each group, each repertoire, each call_type, and each level for unique values.
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
                            assert df[field].sum() == approx(1) or df.empty, f"Expected the sum of"+\
                                f"{field} to be equal to approx 1 or for dataframe to be empty. Instead: \n"+\
                                f"∑({field}) = {df[{field}].sum()} \n"+\
                                f"Variables: \n\tgroup_id:{group_id} \n\trepertoire_id:{rep_id} \n\tcall_type:{combo} \n\tlevel:{level} \n\tmode:{mode} \n\tproductive:{productive}"

def test_field_sum_groups(field):
    '''
    Tests to ensure that the given field variable sums to 1 for each group, each call_type, and each level for unique values.
    The assertion checks if the sum is approximately 1 (using pytest.approx()), or for the dataframe to be empty.
    '''
    for group in groups:
        group_id = group['repertoire_group_id']
        for combo in combos:
            for level in levels:
                filename = f'{group_id}.{processing_stage}.group.{combo}.tsv'
                orig_df = pd.read_csv(os.path.join(curr_dir, filename), delimiter='\t')
                for mode in modes:
                    for productive in productivity:
                        df = orig_df[(orig_df['mode']==mode) & (orig_df['productive']==productive) & (orig_df['level']==level)]
                        assert df[field].sum() == approx(1) or df.empty, f"Expected the sum of"+\
                            f"{field} to be equal to approx 1 or for dataframe to be empty. Instead: \n"+\
                            f"∑({field}) = {df[{field}].sum()} \n"+\
                            f"Variables: \n\tgroup_id:{group_id} \n\tcall_type:{combo} \n\tlevel:{level} \n\tmode:{mode} \n\tproductive:{productive}"

def test_field_sum_repertoires(field):
    '''
    Tests to ensure that the given field variable sums to 1 for each group, each repertoire, each call_type, and each level for unique values.
    The assertion checks if the sum is approximately 1 (using pytest.approx()), or exactly 0.
    '''
    for group in groups:
        group_id = group['repertoire_group_id']
        rep_ids = [rep_obj['repertoire_id'] for rep_obj in group['repertoires']]
        for rep_id in rep_ids:
            for combo in combos:
                for level in levels:
                    filename = f'{rep_id}.{processing_stage}.{combo}.tsv'
                    orig_df = pd.read_csv(os.path.join(curr_dir, filename), delimiter='\t')
                    for mode in modes:
                        for productive in productivity:
                            df = orig_df[(orig_df['mode']==mode) & (orig_df['productive']==productive) & (orig_df['level']==level) & (orig_df['repertoire_id']==rep_id)]
                            assert df[field].sum() == approx(1) or df.empty, f"Expected the sum of"+\
                                f"{field} to be equal to approx 1 or for dataframe to be empty. Instead: \n"+\
                                f"∑({field}) = {df[{field}].sum()} \n"+\
                                f"Variables: \n\tgroup_id:{group_id} \n\trepertoire_id:{rep_id} \n\tcall_type:{combo} \n\tlevel:{level} \n\tmode:{mode} \n\tproductive:{productive}"


# run tests
test_field_sum_repertoires('sequence_frequency')
test_field_sum_repertoires('duplicate_frequency')
test_field_sum_group_repertoires('sequence_frequency')
test_field_sum_group_repertoires('duplicate_frequency')
test_field_sum_groups('sequence_frequency_avg')
test_field_sum_groups('duplicate_frequency_avg')
