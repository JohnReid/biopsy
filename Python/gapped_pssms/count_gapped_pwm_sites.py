#
# Copyright John Reid 2009
#

"""
Code that counts the number of sequences for which a gapped PWM has at least one site in (using varying thresholds).
"""


import logging, sys, pylab as P, numpy as N, hmm, hmm.pssm.logo as L, infpy.roc as roc, cPickle
from optparse import OptionParser
from hmm.pssm import seq_to_numpy, numpy_to_seq
from itertools import imap


def sequences_from_fasta(fasta):
    """Yields sequences from fasta file."""
    import corebio.seq_io.fasta_io
    for seq in corebio.seq_io.fasta_io.iterseq(
        open(fasta, 'r'),
        corebio.seq.reduced_nucleic_alphabet
    ):
        yield seq.description.strip(), (str(seq).strip('nN'), seq.tally())


def load_seqs(filename):
    "Load and convert sequences from fasta file."
    logging.info('Loading sequences: %s', filename)
    sequences = dict(sequences_from_fasta(filename))
    numpy_seqs = dict((desc, hmm.preprocess_sequence(seq_to_numpy(seq))) for desc, (seq, tally) in sequences.iteritems())
    tally = sum(imap(N.array, (tally for desc, (seq, tally) in sequences.iteritems())))
    logging.info('Loaded %d sequences with %d bases', len(sequences), sum(imap(len, (seq for seq, tally in sequences.values()))))
    return numpy_seqs, tally


def build_hmm(freqs, gaps):
    """
    Build a HMM representing the gapped PWM with the given frequencies and gaps. Cannot handle PWMs with consecutive gaps
    or gaps at beginning or end.
    """
    if len(gaps) != len(freqs):
        raise ValueError('Frequencies and gaps must be same length.')
    K = len(gaps)

    # create model
    model = hmm.ModelByStates(M=4, markov_order=0)

    # add background state
    bg = model.add_state()
    bg.pi = model.add_parameter(1.)
    uniform_param = model.add_parameter(.25)
    for m in xrange(bg.M):
        bg.b[m] = uniform_param

    # add the binding site states in positive and negative directions
    positive_states = [model.add_state() for i in xrange(K)]
    negative_states = [model.add_state() for i in xrange(K)]

    # connect background to initial binding site states
    p_binding_site = 0.01
    binding_param = model.add_parameter(p_binding_site/2.)
    not_binding_param = model.add_parameter(1.-p_binding_site)
    bg.add_successor(positive_states[0], binding_param)
    bg.add_successor(negative_states[-1], binding_param)
    bg.add_successor(bg, not_binding_param)
    always_one_param = model.add_parameter(1.)
    positive_states[-1].add_successor(bg, always_one_param)
    negative_states[0].add_successor(bg, always_one_param)

    # set up emissions
    for freq, positive_state, negative_state in zip(freqs, positive_states, negative_states):
        for b, f in enumerate(freq):
            emission_param = model.add_parameter(f)
            positive_state.b[b] = emission_param
            negative_state.b[-b-1] = emission_param
            #positive_state.

    # set up transitions
    for k, gap in enumerate(gaps):
        if gap < 1. and (0 == k or K-1 == k or gaps[k-1] < 1. or gaps[k+1] < 1.):
            raise ValueError('Gaps cannot be at first or last character nor next to another gap.')
        if gap < 1.:
            gap_param = model.add_parameter(gap)
            non_gap_param = model.add_parameter(1.-gap)
            positive_states[k-1].add_successor(positive_states[k], gap_param)
            positive_states[k-1].add_successor(positive_states[k+1], non_gap_param)
            negative_states[k+1].add_successor(negative_states[k-1], non_gap_param)
            negative_states[k+1].add_successor(negative_states[k], gap_param)
        else:
            if 0 != k:
                positive_states[k-1].add_successor(positive_states[k], always_one_param)
            if K-1 != k:
                negative_states[k+1].add_successor(negative_states[k], always_one_param)

    return model


def build_model_by_states(freqs, gaps, p_binding_site=0.001):
    """
    Build a HMM representing the gapped PWM with the given frequencies and gaps. Can handle consecutive gaps
    and gaps at beginning or end.
    """
    if len(gaps) != len(freqs):
        raise ValueError('Frequencies and gaps must be same length.')
    K = len(gaps)

    # normalise frequencies
    freqs = (freqs.T / freqs.sum(axis=1)).T

    # create model
    model = hmm.ModelByStates(M=4, markov_order=0)

    # add background state
    bg = model.add_state()
    bg.pi = model.add_parameter(1.)
    uniform_param = model.add_parameter(.25)
    for m in xrange(bg.M):
        bg.b[m] = uniform_param

    # add the binding site states in positive and negative directions
    positive_states = [model.add_state() for i in xrange(K)]
    negative_states = [model.add_state() for i in xrange(K)]

    # connect background to initial binding site states
    binding_param = model.add_parameter()
    not_binding_param = model.add_parameter(1.-p_binding_site)
    bg.add_successor(positive_states[0], binding_param)
    bg.add_successor(negative_states[-1], binding_param)
    bg.add_successor(bg, not_binding_param)
    always_one_param = model.add_parameter(1.)
    positive_states[-1].add_successor(bg, always_one_param)
    negative_states[0].add_successor(bg, always_one_param)

    # set up emissions
    for freq, positive_state, negative_state in zip(freqs, positive_states, negative_states):
        for b, f in enumerate(freq):
            emission_param = model.add_parameter(f)
            positive_state.b[b] = emission_param
            negative_state.b[-b-1] = emission_param

    # set up transitions
    def setup_transitions(states, gaps):
        for k in xrange(-1, K):
            if -1 == k:
                k_state = bg
                p_skip = p_binding_site/2.
            else:
                k_state = states[k]
                p_skip = 1.
            for m in xrange(k+1, K):
                gap_param = model.add_parameter(p_skip * gaps[m])
                k_state.add_successor(states[m], gap_param)
                p_skip *= (1. - gaps[m])
                if 0. == p_skip:
                    break
            if p_skip > 0.:
                states[k].add_successor(bg, model.add_parameter(p_skip))

    setup_transitions(positive_states, gaps)
    setup_transitions(negative_states[::-1], gaps[::-1])

    return model


def build_hmm_model(freqs, gaps, p_binding_site=.001):
    "@return: A hmm.Model representing the gapped PWM defined by the arguments."
    model_by_states = build_model_by_states(freqs, gaps, p_binding_site=p_binding_site)
    model = hmm.as_model(model_by_states)
    model.normalise()
    return model


def run_on_seqs(model, numpy_seqs):
    """
    Run model on sequences
    """
    total_pos = 0
    total_neg = 0
    num_seqs_with_site = 0
    for desc, seq in numpy_seqs.iteritems():
        LL, states = model.viterbi(seq)
        num_pos = (states == 1).sum()
        num_neg = (states == 1+(model.N-1)/2).sum()
        logging.debug('Ran Viterbi algorithm on %20s: found %3d positive sites and %3d negative sites', desc, num_pos, num_neg)
        total_pos += num_pos
        total_neg += num_neg
        if num_pos + num_neg > 0:
            num_seqs_with_site += 1
    return total_pos, total_neg, num_seqs_with_site



def run_pwm_viterbi(tag, freqs, gaps, positive_seqs, negative_seqs):
    """
    Run the PWM using Viterbi algorithm to classify sequences.
    """
    logging.info('Running PWM: %s', tag)
    logo = L.pssm_as_image(freqs, size=None, transparencies=gaps)
    logo_filename = '%s-logo.png' % tag
    logo.save(logo_filename)
    logging.info('%s: Created logo: %s', tag, logo_filename)
    roc_points = []
    for p_binding in p_binding_params:
        # build model
        model = build_hmm_model(freqs, gaps, p_binding)
        hmm.graph_as_svg(model, '%s-states' % tag, neato_properties={'-Elen':1.4})
        logging.debug('%s: Graphed model', tag)
        pos_total_pos, pos_total_neg, pos_num_seqs_with_site = run_on_seqs(model, positive_seqs)
        logging.debug(
            '%s: p(binding)=%.1e: Positive sequences: Over all sequences: found %4d positive sites and %4d negative sites in %4d/%4d sequences',
            tag,
            p_binding,
            pos_total_pos,
            pos_total_neg,
            pos_num_seqs_with_site,
            len(positive_seqs)
        )
        neg_total_pos, neg_total_neg, neg_num_seqs_with_site = run_on_seqs(model, negative_seqs)
        logging.debug(
            '%s: p(binding)=%.1e: Negative sequences: Over all sequences: found %4d positive sites and %4d negative sites in %4d/%4d sequences',
            tag,
            p_binding,
            neg_total_pos,
            neg_total_neg,
            neg_num_seqs_with_site,
            len(negative_seqs)
        )
        tp = pos_num_seqs_with_site
        fp = neg_num_seqs_with_site
        fn = len(positive_seqs) - pos_num_seqs_with_site
        tn = len(negative_seqs) - neg_num_seqs_with_site
        roc_point = roc.RocCalculator(tp=tp, fp=fp, tn=tn, fn=fn)
        logging.info('%s: p(binding)=%.1e; Specificity=%.3f; Sensitivity=%.3f',
            tag,
            p_binding,
            roc_point.specificity(),
            roc_point.sensitivity(),
        )
        roc_points.append(roc_point)
    return roc_points


def make_classifier(model):
    """
    Given a model, creates a classifier from it. A classifier is a function that is given a sequence and returns the threshold
    above which the sequence would be considered a positive.
    """
    def classifier(sequence):
        "A classifier that takes a sequence and returns at what threshold it would be treated as positive."
        LL, alpha, beta, c = model.forward_backward(sequence)
        alphabeta = alpha * beta
        gamma0 = alphabeta[:,0] / alphabeta.sum(axis=1)
        # return how often we are not in state 0
        return len(gamma0) - gamma0.sum()

    return classifier


def test_hmm_forward_backward(model, seqs):
    """
    Test a HMM on positive and negative sequences using forward-backward algorithm.

    Counts how many bases are expected to be binding sites in each sequence.
    """
    classifier = make_classifier(model)
    scores = map(classifier, seqs)
    scores.sort()
    return scores


def run_pwm_forward_backward(tag, freqs, gaps, positive_seqs, negative_seqs):
    """
    Run the PWM using forward-backward.
    """
    logging.info('Running PWM: %s', tag)
    logo = L.pssm_as_image(freqs, size=None, transparencies=gaps)
    logo_filename = '%s-logo.png' % tag
    logo.save(logo_filename)
    logging.info('%s: Created logo: %s', tag, logo_filename)
    # build model
    model = build_hmm_model(freqs, gaps, .001)
    hmm.graph_as_svg(model, '%s-states' % tag, neato_properties={'-Elen':1.4})
    logging.debug('%s: Graphed model', tag)
    positive_scores = test_hmm_forward_backward(model, positive_seqs.values())
    negative_scores = test_hmm_forward_backward(model, negative_seqs.values())
    return roc.picked_rocs_from_thresholds(positive_scores, negative_scores)


if '__main__' == __name__:
    #
    # Initialise the logging
    #
    logging.basicConfig(level=logging.INFO, format="%(asctime)s-%(name)s:%(levelname)s:%(message)s")
    logging.info(hmm.module_info())



    #
    # Parse the options
    #
    option_parser = OptionParser()
    option_parser.add_option(
        "--min-p-binding-site",
        dest="min_p_binding",
        default=1e-6,
        type='float',
        help="Minimum probability of a binding site in the model."
    )
    option_parser.add_option(
        "--max-p-binding-site",
        dest="max_p_binding",
        default=1e-1,
        type='float',
        help="Maximum probability of a binding site in the model."
    )
    option_parser.add_option(
        "-n",
        "--num-p-binding-params",
        dest="num_p_binding",
        default=10,
        type='int',
        help="Number of p(binding site) parameters to evaluate."
    )
    options, args = option_parser.parse_args()
    for option in option_parser.option_list:
        if option.dest:
            logging.info('%20s: %30s (%s)', option.dest, str(getattr(options, option.dest)), option.help)


    #
    # Load sequences
    #
    if len(args) < 2 or len(args) > 2:
        logging.error('USAGE: %s <positive sequences> <negative sequences>')
    else:
        #
        # Choose p(binding site) params (evenly spaced on a log scale)
        #
        p_binding_params = N.exp(N.linspace(N.log(options.min_p_binding), N.log(options.max_p_binding), options.num_p_binding))

        # import Sp1 matrix definitions
        from sp1 import *


        def test_matrix():
            # TRANSFAC sp1
            N.random.seed(1)
            freqs = N.random.dirichlet(N.ones(4)*.1, size=6)
            gaps = N.array([.1, .5, .4, 1., 1., .3])
            return freqs, gaps

        def test_hmm(tag, pwm):
            freqs, gaps = pwm
            logo = L.pssm_as_image(freqs, size=None, transparencies=gaps)
            logo_filename = '%s-logo.png' % tag
            logo.save(logo_filename)
            logging.info('%s: Created logo: %s', tag, logo_filename)
            model = build_hmm_model(freqs, gaps, .1)
            logging.debug('%s: Created model', tag)
            hmm.graph_as_svg(model, '%s-states' % tag, neato_properties={'-Elen':3.})
            logging.debug('%s: Graphed model', tag)
            return model

        #model = test_hmm('test', test_matrix())
        #glam2_sp1_i7_model = test_hmm('GLAM2-Sp1-i7', glam2_sp1_i7())

        positive_seq_file, negative_seq_file = args
        positive_seqs, tally_positive = load_seqs(positive_seq_file)
        logging.info('Positive tally: %s', str(tally_positive))
        negative_seqs, tally_negative = load_seqs(negative_seq_file)
        logging.info('Negative tally: %s', str(tally_negative))

        def calc_or_unpickle_roc(tag, pwm, positive_seqs, negative_seqs):
            pickle_file = '%s-roc.pickle'%tag
            try:
                roc = cPickle.load(open(pickle_file))
                logging.info('Unpickled ROCs from %s.', pickle_file)
            except:
                logging.info('Could not unpickle %s, calculating from scratch.', pickle_file)
                freqs, gaps = pwm
                roc = run_pwm_forward_backward(tag, freqs, gaps, positive_seqs, negative_seqs)
                cPickle.dump(roc, open(pickle_file, 'wb'))
            return roc

        def do_pwm(tag, pwm, color, marker):
            roc_points = calc_or_unpickle_roc(tag, pwm, positive_seqs, negative_seqs)
            roc.plot_roc_points(roc_points, label=tag, color=color, marker=marker)

        P.figure()
        do_pwm('GLAM2-Sp1-i4', glam2_sp1_i4(), 'magenta', 's')
        do_pwm('GLAM2-Sp1-i7', glam2_sp1_i7(), 'cyan', '^')
        do_pwm('Gapped-Sp1', gapped_sp1(), 'blue', 'o')
        do_pwm('TRANSFAC-Sp1', transfac_sp1(), 'red', 'v')
        do_pwm('MEME-Sp1', meme_sp1(), 'green', 'd')
        roc.plot_random_classifier(label='Random')
        roc.label_plot()
        P.legend(loc='lower right')
        P.savefig('ROC.eps')
        P.savefig('ROC.png')
        P.xlim(0., .2)
        P.ylim(0., .5)
        P.savefig('ROC-zoom.eps')
        P.savefig('ROC-zoom.png')
