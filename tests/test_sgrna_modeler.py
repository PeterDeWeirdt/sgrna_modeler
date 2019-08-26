#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `sgrna_modeler` package."""

import pandas as pd
import os

import numpy as np
from sgrna_modeler import features as ft
from sgrna_modeler import enzymes as en
from sgrna_modeler import activity_data as da

def test_encoding(seqs=None):
    if seqs is None:
        seqs = ['ACT', 'ACT', 'GCT', 'AAT', 'AAA']
    encoded = ft.encode_seqs(seqs)
    np.testing.assert_array_equal(encoded[3, :, :],
                                  np.array([[1, 0, 0, 0],
                                            [1, 0, 0, 0],
                                            [0, 0, 0, 1]]))

def test_context_order():
    assert ft.get_context_order(4) == ['1', '2', '3', '4']

def test_guide_seq():
    # SpCas9
    assert ft.get_guide_sequence('CGTCCCCATCCACGGCCTTCACCCGGGCAG', 5, 20) == 'CCCATCCACGGCCTTCACCC'
    # AsCas12a
    assert ft.get_guide_sequence('CCAGTTTGAACTCTCGCCCATCACCTATCAGTGC', 9, 20) == 'AACTCTCGCCCATCACCTAT'


def test_featurization():
    features = {'Pos. Ind. 1mer',
             'Pos. Ind. 2mer',
             'Pos. Ind. 3mer',
             'Pos. Ind. Zipper',
             'Pos. Dep. 1mer',
             'Pos. Dep. 2mer',
             'Pos. Dep. 3mer',
             'Pos. Dep. Zipper',
             'Pos. Ind. Rep.',
             'GC content',
             'Tm',
             'Cas9 PAM',
             'Physio',
             'OOF Mutation Rate',
             'Double Zipper'}

    kmers = pd.Series(['ACTGGTGGG'])
    one_hot = ft.featurize_guides(kmers, features, guide_start=3, guide_length=6).iloc[0]
    print(one_hot)
    assert (one_hot['3T'] == 1)
    assert (one_hot['TGG'] == 0.5)
    assert (one_hot['GC content'] == 2/3)
    assert (one_hot['bendability'] != 0)
    assert (one_hot['Tm, context'] != 0)
    assert (one_hot['1AC'] == 1)
    assert (one_hot['7GGG'] == 1)
    assert (one_hot['G'] == 2/3)
    assert (one_hot['TG'] == 0.4)
    assert (one_hot['1T'] == 0)

def test_pam_splits():
    split1 = en.get_pam_splits(['TTTV'])
    assert set(split1) == {'TTTA', 'TTTC', 'TTTG'}
    split2 = en.get_pam_splits(['RGG', 'SGG'])
    assert (split2) == {'AGG', 'GGG', 'CGG'}

def test_data_load():
    fname = os.path.join(os.path.dirname(__file__), 'kim_2018_test.csv')
    data = pd.read_csv(fname)
    cas12a = en.Enzyme(guide_start=9, guide_length=23,
                       pam_start=5, pams = ['TTTN'], context_length=34)
    kim_2018_train = da.Activity_Data(data, cas12a, 'Context Sequence',
                                      'Indel frequency')
    x, y = kim_2018_train.get_xy()
    assert x.shape == y.shape
    assert kim_2018_train.enzyme.context_length == 34
