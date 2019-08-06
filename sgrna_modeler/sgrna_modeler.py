# -*- coding: utf-8 -*-
import numpy as np
"""Main module."""
nt_codes = {'A':[1,0,0,0],
            'C':[0,1,0,0],
            'G':[0,0,1,0],
            'T':[0,0,0,1]}
def encode_seqs(seqs):
    # 3d array with samples x position x nt
    encoded_seqs = np.array([[nt_codes.get(x) for x in seq] for seq in seqs])
    return encoded_seqs

def get_nt_count(encoded_seqs, nt):
    counts = np.sum(encoded_seqs[:,:,nt_codes.get(nt).index(1)], axis = 1)
    return counts

