#
# Copyright John Reid 2006-2009
#

"""
Code to deal with binding site test cases.
"""


from _biopsy import *
from binding_hit import *

import random

def calc_sensitivity( tp, fn ):
    return tp + fn and float(tp)/(float(tp)+float(fn)) or 0

def calc_specificity( tn, fp ):
    return tn + fp and float(tn)/(float(tn)+float(fp)) or 0

def evaluate_discriminator_at_threshold(
        discriminators,
        truths,
        threshold
):
    """Evaluates a discriminator at a given threshold and returns
    (sensitivity, specificity)"""
    tp = fp = fn = tn = 0
    for d, should_be_positive in zip(discriminators, truths):
        is_positive = d >= threshold
        if should_be_positive and is_positive: tp += 1
        elif should_be_positive and not is_positive: fn += 1
        elif not should_be_positive and is_positive: fp += 1
        elif not should_be_positive and not is_positive: tn += 1
    sensitivity = calc_sensitivity(tp,fn)
    specificity = calc_specificity(tn,fp)
    return ( sensitivity, specificity )

def add_negatives_to_test_case( tc, pssm_names ):
    """Generates random negative test cases"""
    # hash the positives to generate seed for random
    rand = random.Random( "".join( tc.positives ).__hash__() )
    tc.negatives = SequenceVec()
    while len(tc.negatives) != len(tc.positives):
        p = rand.choice(pssm_names)
        if p not in tc.positives:
            tc.negatives.append(p)

TestCase.add_negatives_from = add_negatives_to_test_case

#
# Our tests are passed around as tuples:
# test[0] is the test case
# test[1] is whether this is a positive test
#

class Test:
    def __init__( self, tc, truth ):
        self.tc = tc
        self._truth = truth
    def truth( self ):
        return self._truth
    def pssms( self ):
        return self.truth() and self.tc.positives or self.tc.negatives
    def sequences( self ):
        return self.tc.sequences

def test_case_truth( test ):
    """Should a pssm be a positive or a negative test case input?"""
    return test.truth()

class PhyloClassifier:
    """Classifies tests. Returns a prob [0,1] for whether the pssm is predicted
    to be in the test case"""
    def __init__( self, threshold, phylo_threshold ):
        self._threshold = threshold
        self._phylo_threshold = phylo_threshold
        self._p_bindings = { }
    def __call__( self, test ):
        hits = score_pssms_on_phylo_sequences(
                test.pssms(),
                test.sequences(),
                self._threshold,
                self._phylo_threshold )[ 0 ]
        p_bindings = get_max_p_binding_over_hits( hits )
        return (
                0 != len(p_bindings)
                and max(p_bindings.values())
                or 0.0
        )

class NonPhyloClassifier:
    """Classifies tests. Returns a prob [0,1] for whether the pssm is predicted
    to be in the test case"""
    def __init__( self, threshold ):
        self._threshold = threshold
        self._p_bindings = { }
    def __call__( self, test ):
        hits = score_pssms_on_sequence(
                test.pssms(),
                test.sequences()[0],
                self._threshold )
        p_bindings = get_max_p_binding_over_hits( hits )
        result = (
                0 != len(p_bindings)
                and max(p_bindings.values())
                or 0.0
        )
        # print test.truth(), result, test.pssms()
        return result

class TransfacClassifier:
    """Classifies tests. Returns a prob [0,1] for whether any pssm is predicted
    to be in the test case"""
    def __init__( self, threshold_fn = get_transfac_pssm_min_fp_threshold ):
        self._threshold_fn = threshold_fn
    def __call__( self, test ):
        for pssm in test.pssms():
            p = get_pssm(pssm).pssm
            l = len(p)
            try:
                threshold = self._threshold_fn( pssm )
            except:
                print '********** Cannot get threshold for:', pssm, '**********'
                continue
            for i in range( len( test.sequences()[0] ) - l + 1 ):
                s = test.sequences()[0][i:i+l]
                if threshold <= score_pssm( p, s ):
                    return 1.0
                s_rc = reverse_complement( s )
                if threshold <= score_pssm( p, s_rc ):
                    return 1.0
        # print s, test_case.sequences[0][-l:]
        return 0.0

def generate_test_inputs( test_cases, max_num_sequences = 6 ):
    """A generator for test inputs"""
    # for each test case
    for tc in test_cases:
        print '# seqs: %d, # pssms: %d, seq lengths: %s' % (
                len( tc.sequences ),
                len( tc.positives ),
                ','.join( str( len( s ) ) for s in tc.sequences )
        )

        if len( tc.sequences ) > max_num_sequences:
            print 'Too many sequences, skipping...'
            continue

        yield Test( tc, True )
        yield Test( tc, False )


def evaluate_test_case_classifier(
        test_cases,
        classifier,
        truth,
        num_thresholds = 100
):
    """Analyses the test cases"""

    discriminators = [ classifier( i ) for i in generate_test_inputs( test_cases ) ]
    truths = [ truth( i ) for i in generate_test_inputs( test_cases ) ]
    print '# tests:', len(truths)
    # for i in range(len(truths)):
    #       print truths[i], discriminators[i]

    # for each threshold
    roc_points = { }
    for t in [ float(i) / ( num_thresholds - 1 ) for i in range(num_thresholds) ]:

        # Sensitivity and specificity
        (sensitivity, specificity) = \
                evaluate_discriminator_at_threshold(
                        discriminators,
                        truths,
                        t )

        # print sensitivity, specificity
        roc_points[ t ] = (
                sensitivity,
                1.0 - specificity
        )

    return roc_points

def plot_roc_points(
        roc_points,
        **keywords
):
    from pylab import plot
    x = []
    y = []
    thresholds = roc_points.keys()
    thresholds.sort()
    for t in thresholds:
        # print t, point
        y.append( roc_points[t][0] )
        x.append( roc_points[t][1] )
    # print x
    # print y
    x.sort()
    y.sort()
    plot( x, y, **keywords )
