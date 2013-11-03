#
# Copyright John Reid 2010, 2011
#

import _biopsy as biopsy
import _bifa as bifa
import numpy as N

print 'BiFa C++ module version %s' % bifa.version()

# make a PSSM
pssm = N.log(
    N.array((
        (.7, .1, .1, .1),
        (.1, .7, .1, .1),
        (.1, .1, .7, .1),
        (.1, .1, .1, .7),
        (.25, .25, .25, .25)
    ))
)
print pssm

# make sure score is correct when given int sequence
assert 3.05 == bifa.score_word(N.exp(pssm), [0, 1, 2, 3, 4])

# make sure score is correct when given string sequence
assert 3.05 == bifa.score_word(N.exp(pssm), 'acgtn')

# try too short sequence
try: bifa.score_word(pssm, [1, 1, 1, 1])
except RuntimeError: pass
else: assert not "Should have thrown error"

# hide non-int in sequence
try: bifa.score_word(pssm, [1, .25, 1, 1, 1])
except TypeError: pass
else: assert not "Should have thrown error"

# try non-int sequence
try: bifa.score_word(pssm, [.25, .25, .25, .25, .25])
except RuntimeError: pass
else: assert not "Should have thrown error"


print bifa.score_one_strand(pssm, [0, 1, 2, 3, 4, 0, 1, 2, 3, 4])
for hit in bifa.score_sequence(pssm, [0, 1, 2, 3, 4, 0, 1, 2, 3, 4]):
    print hit.p_binding, hit.location.position, hit.location.positive_strand
