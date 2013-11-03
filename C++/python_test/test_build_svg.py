#
# Copyright John Reid 2010
#

import _biopsy as B
B.init()

B.PssmParameters.use_score = True

seq = "ACGCGAGCAGGGTCATTAAATCGAGCGTCGCGGCGCGCGACAAGGACGGCATTATTAGCGTGCTACGACTACGACTTG"
pssms = B.SequenceVec()
pssms.extend(('M00023', 'M00436', 'M00716', 'M00938', 'R04653'))
hits = B.score_pssms_on_sequence(pssms, seq)

svg_args = B.BuildSvgArgs(file='test.svg')
B.build_analysis_svg(hits, None, seq, svg_args)
