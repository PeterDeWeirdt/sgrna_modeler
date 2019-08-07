#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `sgrna_modeler` package."""

import pytest

import numpy as np
from sgrna_modeler import sgrna_modeler as sm

seqs = ['ACT','ACT','GCT','AAT','AAA']

def test_encoding():
    encoded = sm.encode_seqs(seqs)
    np.testing.assert_array_equal(encoded[3,:,:],
                                  np.array([[1,0,0,0],
                                            [1,0,0,0],
                                            [0,0,0,1]]))

def test_nt_counts():
    encoded = sm.encode_seqs(seqs)
    A_content = sm.get_nt_count(encoded, 'A')
    np.testing.assert_array_equal(A_content, np.array([1,1,0,2,3]))


