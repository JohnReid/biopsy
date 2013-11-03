#
# Copyright John Reid 2010
#

import _biopsy as B
B.init()

seq = "ACGCGAGCAGCATCATTATATCGAGCGACGCGGCGCGCGACAAGGACGGCATTATTAGCGAGCTACGACTACGACTTG"
pssms = B.SequenceVec()
hits = B.HitVec()
p_binds = B.score_pssm_on_sequence('M00023', seq, .03, hits)
