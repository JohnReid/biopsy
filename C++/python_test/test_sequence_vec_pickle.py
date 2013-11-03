import _biopsy as B
import cPickle

sequences = B.SequenceVec()
sequences.append('acgt')
sequences.append('tttt')
sequences.append('gggg')

sequences_copy = cPickle.loads(cPickle.dumps(sequences))

for s, c in zip(sequences, sequences_copy):
    assert s == c
