# Enzymes

nt_codes = {'A': 'A',
            'C': 'C',
            'T': 'T',
            'G': 'G',
            'R': ['A', 'G'],
            'Y': ['C', 'T'],
            'S': ['G', 'C'],
            'M': ['A', 'C'],
            'K': ['G', 'T'],
            'B': ['G', 'C', 'T'],
            'H': ['A', 'C', 'T'],
            'D': ['A', 'G', 'T'],
            'V': ['A', 'G', 'C'],
            'N': ['A', 'C', 'T', 'G']}

def expand_seq(encoded_seq, seq_str='', seqs=[]):
    if not encoded_seq:
        seqs.append(seq_str)
        return seqs
    else:
        nts = nt_codes[encoded_seq[0]]
        for nt in nts:
            base = seq_str
            base += nt
            expand_seq(encoded_seq[1:], base, seqs)
    return seqs

def get_pam_splits(pam_list):
    #   Take a list of PAM string and return a list of PAMs
    all_pams = set()
    for pam in pam_list:
        all_pams |= set(expand_seq(pam, seq_str='', seqs=[]))
    return all_pams

class Enzyme(object):
    # Default paramaters are for AsCas12a
    def __init__(self, guide_start = 9, guide_length = 23,
                 pam_start = 5, pams = ['TTTN'], context_length = 34):
        self.guide_start = guide_start
        self.guide_length = guide_length
        self.pam_start = pam_start
        self.pams = pams
        self.context_length = context_length

cas9 = Enzyme(context_length=30, guide_length=20, guide_start=5,
              pam_start=25, pams=['TTTV'])
