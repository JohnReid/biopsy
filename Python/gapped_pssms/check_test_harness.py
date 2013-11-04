#
# Copyright John Reid 2009
#

"""
Code to test how test harness reacts to sequences with Ns
"""

import cPickle, os, sys, logging, numpy, corebio.seq, hmm, numpy as N
from itertools import chain
from test_harness_2 import make_classifier
from parse_gapped_format import parse_models, build_hmm_from_semi_parsed

old_results_dir = '/home/reid/Analysis/GappedPssms/test-harness/Results/2009-07-30'
new_results_dir = '/home/reid/Analysis/GappedPssms/test-harness/Results/2009-08-01'

#meme_neg=list(chain(*[cPickle.load(open('T00671-x%d-neg-r1-back-MEME.results'%i)) for i in xrange(1,6)]))
#meme_neg.sort()
#meme_pos = list(chain(*[cPickle.load(open('T00671-x%d-MEME.results'%i)) for i in xrange(1,6)]))
#meme_pos.sort()

old_seq="""tacatggtacttgctttaaaaagtagaatagggagtatgaggttactttatgtgaagatactcccgctgaacaataaata
aatactattactatacccacgacctccagaaattcactggataaccagtaagacaacttctactcatttcttcatattcc
tacttattcaagttgtagccttcatagttgataaaaaatcagcacacattaagaaaacaataacagaactattttcttca
catgacttttattccttaatccagactgttaaaaggactgcaagacaaattgtttttcaatcagatttttttctccacca
gatgtctatgtgaatttcatattgttttagacaaaaatgctcattccttcggtctaagtactatgtcatattttgttttt
tcaagccttcaaattttgtgctggtggttacttcatatacattctatggttaatctttaaagagaagttttaaaagtctg
attcaaaatttcagttcactcgctatgtattttaaaaattaaaatttatgaaattcaattttaaaaatctaaaagttatc
taaaaaggtctatgacttatcaaatttcaataagctgactgttagcagtattaaaaaatattaaatatgctaacagtaaa
aatcatgaatacacattaggcatttaatatgtatctggcaaatttgaatacataaagggaataggcagagttcacagatt
aatatttcttacctctacaataagaagaaataccttgttctatgagcagctgccatactttcagacatgtttctgacttt
tagataattaacaaatcctctgaagaaaaggagcaggcctgagaaggttgaaataatatggatatactatgtttttatac
agaaaagggcaagataaatttaaagtagacaattataaacaagatcctcttggcttaagga""".replace('\n', '')

new_seq="""TACATGGTACTTGCTTTAAAAAGTAGAATAGGGAGTATGAGGTTACTTTATGTGAAGATACTCCCGCTGAACAATAAATA
AATACTATTACTATACCCACGACCTCCAGAAATTCACTGGATAACCAGTAAGACAACTTCTACTCATTTCTTCATATTCC
TACTTATTCAAGTTGTAGCCTTCATAGTTGATAAAAAATCAGCACACATTAAGAAAACAATAACAGAACTATTTTCTTCA
CATGACTTTTATTCCTTAATCCAGACTGTTAAAAGGACTGCAAGACAAATTGTTTTTCAATCAGATTTTTTTCTCCACCA
GATGTCTATGTGAATTTCATATTGTTTTAGACAAAAATGCTCATTCCTTCGGTCTAAGTACTATGTCATATTTTGTTTTT
TCAAGCCTTCAAATTTTGTGCTGGTGGTTACTTCATATACATTCTATGGTTAATCTTTAAAGAGAAGTTTTAAAAGTCTG
ATTCAAAATTTCAGTTCACTCGCTATGTATTTTAAAAATTAAAATTTATGAAATTCAATTTTAAAAATCTAAAAGTTATC
TAAAAAGGTCTATGACTTATCAAATTTCAATAAGCTGACTGTTAGCAGTATTAAAAAATATTAAATATGCTAACANNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNATACATAAAGGGAATAGGCAGAGTTCACAGATT
AATATTTCTTACCTCTACAATAAGAAGAAATACCTTGTTCTATGAGCAGCTGCCATACTTTCAGACATGTTTCTGACTTT
TAGATAATTAACAAATCCTCTGAAGAAAAGGAGCAGGCCTGAGAAGGTTGAAATAATATGGATATACTATGTTTTTATAC
AGAAAAGGGCAAGATAAATTTAAAGTAGACAATTATAAACANNNNNNNNNNNNNNNNNGGA""".replace('\n', '')

def convert_seq(seq):
    return numpy.array(corebio.seq.Seq(seq, alphabet=corebio.seq.reduced_nucleic_alphabet).ords())

old_pp = hmm.preprocess_sequence(convert_seq(old_seq))
new_pp = hmm.preprocess_sequence(convert_seq(new_seq))

#meme_dir = '/home/reid/Analysis/GappedPssms/MEME/x-validate'
#pssm_file = os.path.join(meme_dir, 'T00671-1.pssm')
pssm_file = '/home/john/Analysis/GappedPssms/MEME/x-validate/vm-T00671-motif-h2-v9-x1.pssm'
semi_parsed_models = list(parse_models(open(pssm_file)))
if len(semi_parsed_models) > 1:
    print >> sys.stderr, 'For the moment we can only handle one model at a time.'
    sys.exit(-1)
parsed = semi_parsed_models[0]
logging.info(str(parsed))
model, traits = build_hmm_from_semi_parsed(parsed)
classifier = make_classifier(model)

def test_seq(seq):
    return classifier(convert_seq(seq))

print 'Old sequence (without Ns):', classifier(old_pp)
print 'New sequence (with    Ns):', classifier(new_pp)

LL, alpha, beta, c = model.forward_backward(new_pp)
alphabeta = alpha * beta
gamma0 = alphabeta[:,0] / alphabeta.sum(axis=1)
gamma0[N.where(new_pp.as_numpy()==4)] = 1.
