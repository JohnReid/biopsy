#
# Copyright John Reid 2009
#

"""
Code that counts the number of sequences for which a gapped PWM has at least one site in (using varying thresholds).
"""


import logging, pylab as P, numpy as N, hmm, hmm.pssm.logo as L, infpy.roc as roc, cPickle, infpy.roc as R, os
from optparse import OptionParser
from sp1 import *
from count_gapped_pwm_sites import test_hmm_forward_backward, sequences_from_fasta, seq_to_numpy, build_hmm_model
from itertools import imap
from cookbook.dicts import ChoiceDict

def load_seqs(filename):
    "Load and convert sequences from fasta file."
    logging.info('Loading sequences: %s', filename)
    sequences = dict(sequences_from_fasta(filename))
    numpy_seqs = dict((desc, hmm.preprocess_sequence(seq_to_numpy(seq))) for desc, (seq, tally) in sequences.iteritems())
    tally = sum(imap(N.array, (tally for desc, (seq, tally) in sequences.iteritems())))
    logging.info('Loaded %d sequences with %d bases', len(sequences), sum(imap(len, (seq for seq, tally in sequences.values()))))
    logging.info(str(tally))
    return numpy_seqs


#
# Initialise the logging
#
logging.basicConfig(level=logging.INFO, format="%(asctime)s-%(name)s:%(levelname)s:%(message)s")
logging.info(hmm.module_info())


#
# Parse options
#
option_parser = OptionParser()
option_parser.add_option(
    '--num-roc-points',
    dest='num_points',
    type='int',
    default=24,
    help="Number of points on each ROC curve."
)
options, args = option_parser.parse_args()
for option in option_parser.option_list:
    if option.dest:
        logging.info('%30s: %30s (%s)', option.dest, str(getattr(options, option.dest)), option.help)


class Sequences(dict):
    def __missing__(self, key):
        value = load_seqs(sequence_filenames[key])
        self[key] = value
        return value

sequence_filenames = {
    'positive' : '/home/reid/Data/GappedPssms/Full-Sp1/Full-Sp1.fasta',
    'shuffled' : '/home/reid/Data/GappedPssms/Full-Sp1/Full-Sp1-shuffled.fasta',
    'r1-back'  : '/home/reid/Data/GappedPssms/Full-Sp1/vm-r1-back-human-hsNCBI36-v4-big-shortened.fa',
    'r3-TSS'   : '/home/reid/Data/GappedPssms/Full-Sp1/vm-r3-TSS-human-hsNCBI36-v4-big-shortened.fa',
}
#sequence_filenames = {
#    'positive' : '/home/reid/Data/GappedPssms/Full-Sp1/Test-Sp1-shuffled.fasta',
#    'test'     : '/home/reid/Data/GappedPssms/Full-Sp1/Test-Sp1.fasta',
#}

sequences = Sequences()
backgrounds = set(sequence_filenames.keys())
backgrounds.remove('positive')
scores = dict()
methods = args
sp1_pssms = all_sp1_pssms()
for tag in methods:
    score_pickle_file = '%s-scores.pickle' % tag
    try:
        positive_scores, negative_scores = cPickle.load(open(score_pickle_file))
        logging.info('%s: Unpickled ROCs from %s.', tag, score_pickle_file)
    except:
        logging.info('%s: Could not ROCs from unpickle %s, calculating from scratch.', tag, score_pickle_file)
        freqs, gaps = sp1_pssms[tag]
        freqs = (freqs.T / freqs.sum(axis=1)).T
        logo = L.pssm_as_image(freqs, size=None, transparencies=gaps)
        logo_filename = '%s-logo.png' % tag
        logo.save(logo_filename)
        logging.info('%s: Created logo: %s', tag, logo_filename)
        model = build_hmm_model(freqs, gaps, .001)
        hmm.graph_as_svg(model, '%s-states' % tag, neato_properties={'-Elen':1.4})
        logging.debug('%s: Graphed model', tag)
        positive_scores = test_hmm_forward_backward(model, sequences['positive'].values())
        negative_scores = dict(
            (bg, test_hmm_forward_backward(model, sequences[bg].values()))
            for bg in backgrounds
        )
        cPickle.dump((positive_scores, negative_scores), open(score_pickle_file, 'wb'))
    scores[(tag,)] = positive_scores
    for bg, score in negative_scores.iteritems():
        scores[(tag, bg)] = score


#
# Generate ROCs
#
markers = ChoiceDict(('s', '^', 'v', 'o', 'h', 'd', '+', 'x'))
def name(method):
    if 'Gapped' == method:
        return 'Gapped-old'
    if 'Gapped-new' == method:
        return 'Gapped'
    if 'Ungapped-new' == method:
        return 'Ungapped'
    if 'GLAM2-i7' == method:
        return 'GLAM2'
    return method

for bg in backgrounds:

    # ROC curves
    P.figure()
    for method in methods:
        rocs = R.picked_rocs_from_thresholds(
            scores[(method,)],
            scores[(method, bg)],
            num_points=options.num_points
        )
        auc = R.area_under_curve(rocs)
        R.plot_roc_points(rocs, label='%.2f %s'%(auc,name(method)), marker=markers[method])
    R.plot_random_classifier(label='0.50 Random')
    R.label_plot()
    P.legend(loc='lower right')
    P.title('Full Sp1 - %s' % bg)
    P.savefig('ROC-Sp1-%s.png' % bg)
    P.savefig('ROC-Sp1-%s.eps' % bg)

    # precision-recall curves
    P.figure()
    for method in methods:
        rocs = R.picked_rocs_from_thresholds(
            scores[(method,)],
            scores[(method, bg)],
            num_points=options.num_points
        )
        R.plot_precision_versus_recall(rocs, label=name(method), marker=markers[method])
    R.label_precision_versus_recall()
    P.legend(loc='lower left')
    P.title('Full Sp1 - %s' % bg)
    P.savefig('Precision-Recall-Sp1-%s.png' % bg)
    P.savefig('Precision-Recall-Sp1-%s.eps' % bg)
