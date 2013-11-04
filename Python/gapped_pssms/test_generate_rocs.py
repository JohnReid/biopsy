#
# Copyright John Reid 2009
#

"""
Code to generate ROC curves. Designed to work with test_harness_2.
"""

import logging, os, pylab as P, infpy.roc as R, numpy as N
from test_harness_2 import TestHarness, fragment_name
from optparse import OptionParser
from cookbook.dicts import DictOf


def reduce_sequence(sequence, desired_length):
    """Reduces a sequence to the desired length by removing some of its elements uniformly."""
    if len(sequence) < desired_length:
        raise RuntimeError('Cannot reduce sequence to longer length.')
    indexes = N.arange(desired_length) * len(sequence) / desired_length
    return [sequence[i] for i in indexes]



#
# Initialise the logging
#
logging.basicConfig(level=logging.INFO)


#
# Parse options
#
option_parser = OptionParser()
TestHarness.add_options(option_parser)
option_parser.add_option(
    '--num-roc-points',
    dest='num_points',
    type='int',
    default=24,
    help="Number of points on each ROC curve."
)
option_parser.add_option(
    '--num-negative',
    dest='num_negative',
    type='int',
    default=50,
    help="Number of negative examples used for AUC50 calculation."
)
options, args = option_parser.parse_args()
for option in option_parser.option_list:
    if option.dest:
        logging.info('%30s: %30s (%s)', option.dest, str(getattr(options, option.dest)), option.help)


#
# Parse arguments
#
methods = args
logging.info('Methods are: %s', ' '.join(methods))
if not len(methods):
    raise RuntimeError('No methods specified on command line')


#
# Set up the test harness
#
harness = TestHarness(options)



#
# Merge scores
#
scores = DictOf(list) # indexed by (method, fragment) or (method, fragment, bg)
def add_results(key, results):
    scores[key] += results
    scores[key].sort()
for method in methods:
    for fragment in harness.options.fragments:
        for fold in harness.folds():
            dataset = (fragment, fold)
            results = harness.results(dataset, method)
            add_results((method, fragment), results)
            add_results((method, 'Overall'), results)
            for bg in harness.options.backgrounds:
                dataset = (fragment, fold, bg)
                results = harness.results(dataset, method)
                add_results((method, fragment, bg), results)
                add_results((method, 'Overall', bg), results)


#
# Calculate 'Overall-normed' results in which each fragment has same weight. I.e. same number of sequences
#
smallest_dataset = min(len(scores[(methods[0], fragment)]) for fragment in harness.options.fragments)
logging.info('Smallest dataset has length %d', smallest_dataset)
for method in methods:
    for fragment in harness.options.fragments:
        add_results((method, 'Overall-normalised'), reduce_sequence(scores[(method, fragment)], smallest_dataset))
        for bg in harness.options.backgrounds:
            add_results((method, 'Overall-normalised', bg), reduce_sequence(scores[(method, fragment, bg)], smallest_dataset))


class ChoiceDict(dict):
    """
    Dictionary that assigns one value in a list to each key.
    """

    def __init__(self, values, cycle=False):
        "Construct with given values."
        dict.__init__(self)
        import itertools
        if cycle:
            self.it = itertools.cycle(values)
            "Iterator over values."
        else:
            self.it = iter(values)
            "Iterator over values."

    def __missing__(self, key):
        "Called when key is missing."
        self[key] = self.it.next()
        return self[key]


markers = ChoiceDict(('s', '^', 'v', 'o', '+', 'd', 'h', 'x'))
colors = ChoiceDict(('blue', 'green', 'red', 'cyan'))
#
# Generate ROCs
#
fragments = ['Overall', 'Overall-normalised']
fragments += harness.options.fragments
num_fragments = len(harness.options.fragments)
for fragment in fragments:
    aucs = dict()
    for bg in harness.options.backgrounds:
        num_negative = fragment.startswith('Overall') and options.num_negative * num_fragments or options.num_negative # use more negative examples in overall results
        roc_thresholds = dict(
            (
                method,
                R.create_rocs_from_thresholds(
                    scores[(method, fragment)],
                    scores[(method, fragment, bg)],
                    num_points=options.num_points
                )
            ) for method in methods
        )
        auc50s = dict(
            (
                method,
                R.auc50(
                    scores[(method, fragment)],
                    scores[(method, fragment, bg)],
                    num_negative=num_negative,
                    num_points=options.num_points
                )
            ) for method in methods
        )
        aucs[bg] = dict()
        aucs[bg]['AUC'] = dict(
            (method, R.area_under_curve([roc for roc, t in roc_thresholds[method]]))
            for method in methods
        )
        aucs[bg]['AUC50'] = dict(
            (method, auc50s[method][0])
            for method in methods
        )


        # ROC curves
        P.figure()
        for method in methods:
            rocs = [roc for roc, t in roc_thresholds[method]]
            auc = aucs[bg]['AUC'][method]
            auc50 = aucs[bg]['AUC50'][method]
            R.plot_roc_points(rocs, label=method, marker=markers[method], color=colors[method])
        R.plot_random_classifier(label='Random')
        R.label_plot()
        P.legend(loc='lower right')
        P.title('%s - %s' % (fragment_name(fragment), bg))
        P.savefig(os.path.join(options.results_dir, 'ROC-%s-%s.png' % (fragment, bg)))
        P.savefig(os.path.join(options.results_dir, 'ROC-%s-%s.eps' % (fragment, bg)))
        P.close()


        # precision-recall curves
        P.figure()
        for method in methods:
            rocs = [roc for roc, t in roc_thresholds[method]]
            R.plot_precision_versus_recall(rocs, label=method, marker=markers[method], color=colors[method])
        R.label_precision_versus_recall()
        P.legend(loc='lower left')
        P.title('%s - %s' % (fragment_name(fragment), bg))
        P.savefig(os.path.join(options.results_dir, 'Precision-Recall-%s-%s.png' % (fragment, bg)))
        P.savefig(os.path.join(options.results_dir, 'Precision-Recall-%s-%s.eps' % (fragment, bg)))
        P.close()

    # do AUC bar-chart
    #P.rcParams['xtick.direction'] = 'out'
    def x(a, b, m):
        return b*(len(methods)+1) + m
    P.figure(figsize=(14,6))
    for a, auc in enumerate(('AUC', 'AUC50')):
        ax = P.subplot(1, 2, a+1)
        xlocs = []
        xlabels = []
        for b, bg in enumerate(harness.options.backgrounds):
            rects = P.bar(
                [x(a, b, m) for m in xrange(len(methods))],
                [aucs[bg][auc][method] for method in methods],
                color=[colors[method] for method in methods],
                width=1.
            )
            xlocs.append(x(a, b, len(methods)/2.))
            xlabels.append(bg)
        P.xticks(xlocs, xlabels)
        P.ylim(0,1)
        P.title(auc)
        for tl in ax.get_xticklines():
            tl.set_visible(False)
    #P.text(.5, .9, '%s' % fragment, horizontalalignment='center')
    P.figlegend(rects, methods, loc='upper right')
    P.savefig(os.path.join(options.results_dir, 'AUC-%s.png' % fragment))
    P.savefig(os.path.join(options.results_dir, 'AUC-%s.eps' % fragment))
    P.close()
