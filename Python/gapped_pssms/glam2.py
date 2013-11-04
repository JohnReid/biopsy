#
# Copyright John Reid 2009
#

"""
Code to handle GLAM2 output files.
"""

import numpy as N, glob, re, os, hmm.pssm.logo as L, logging
from itertools import imap, repeat
from cookbook.dicts import DictOf

_logger = logging.getLogger(__name__)
#_logger.setLevel(logging.INFO)


class GLAM2Output(object):
    """
    Collates information from GLAM2 output file.
    """

    def __str__(self):
        "String representation."
        return """
Sequences: %d
Greatest sequence length: %d
Residue counts: %s

Score: %.2f  Columns: %d  Sequences: %d

   a    c    g    t   Del Ins Score
%s
""" % (
        self.input_sequences,
        self.biggest_length,
        ' '.join('%s=%d' % (residue, count) for residue, count in self.residue_counts.iteritems()),
        self.score, self.columns, self.hit_sequences,
        '\n'.join(
            '%s %5d     %f' % (
                ' '.join('%4d' % int(f) for f in freq),
                gap,
                score
            ) for freq, gap, score in zip(self.freqs, self.gaps, self.pos_scores)
        )
    )


    def __cmp__(self, other):
        "Compare to other GLAM2Output. One with higher score is greater."
        return cmp(self.score, other.score)


    @staticmethod
    def parse(f):
        "Parse GLAM2 output."
        result = GLAM2Output()
        for l in f:
            l = l.strip()
            if l.startswith('Sequences:'):
                result.input_sequences = int(l.split(':')[1])
            elif l.startswith('Greatest sequence length:'):
                result.biggest_length = int(l.split(':')[1])
            elif l.startswith('Residue counts:'):
                fields = l.split(' ')
                def convert_field(field):
                    res, count = field.split('=')
                    return res, int(count)
                result.residue_counts = dict(map(convert_field, fields[2:]))
            elif l.startswith('Score:'):
                fields = l.split()[1::2]
                result.score = float(fields[0])
                result.columns = int(fields[1])
                result.hit_sequences = int(fields[2])
            elif l.startswith('a  c  g  t Del Ins Score'):
                result.freqs = []
                result.gaps = []
                result.pos_scores = []
                for l in f:
                    fields = l.split()
                    if 6 == len(fields):
                        result.freqs.append(map(int, fields[:4]))
                        result.gaps.append(int(fields[4]))
                        result.pos_scores.append(float(fields[5]))

        for attr in [
            'score',
            'freqs',
            'gaps',
        ]:
            if not hasattr(result, attr):
                raise RuntimeError('Parsing unsuccessful - no "%s" attr' % attr)

        return result

    def freqs_and_gaps(self):
        "@return: (freqs, gaps)"
        freqs = N.array(self.freqs, dtype=float)
        freqs = (freqs.T / freqs.sum(axis=1)).T
        gaps = N.array(self.gaps, dtype=float) / self.hit_sequences
        return freqs, 1.-gaps


def pick_best(filenames):
    "Pick the best output from the output files given."
    outputs = [(output, filename) for output, filename in zip(map(GLAM2Output.parse, imap(open, filenames)), filenames)]
    outputs.sort()
    return outputs[-1][1]


def make_tag(fragment, cross_fold_index, seed):
    "@return: A tag representing the arguments."
    return '%s-x%d-s%d' % (fragment, cross_fold_index, seed)


def output_filename(fragment, cross_fold_index, seed):
    "@return: The filename of the output file for these arguments."
    return 'glam2-%s.out' % make_tag(fragment, cross_fold_index, seed)


_output_file_re = re.compile('glam2-(.*)-x([0-9]+)-s([0-9]+).out')


def interpret_output_filename(filename):
    "Parse the information in the filename."
    match = _output_file_re.search(filename)
    fragment, cross_fold_index, seed = match.groups()
    cross_fold_index = int(cross_fold_index)
    seed = int(seed)
    return fragment, cross_fold_index, seed


def group_output_files(directory):
    "Find all glam2 output files in a directory and group them by run."
    all_files = glob.glob(os.path.join(directory, 'glam2*.out'))
    basenames = map(os.path.basename, all_files)
    #logging.info('\n'.join(basenames))
    result = DictOf(set)
    for filename in basenames:
        fragment, cross_fold_index, seed = interpret_output_filename(filename)
        #logging.info(fragment, cross_fold_index, seed)
        result[(fragment, cross_fold_index)].add(seed)
    return result


def best_files(directory):
    "Work out which seeds have given best motifs out of all output files in given directory."
    output_groups = group_output_files(directory)
    for (fragment, cross_fold_index), seeds in output_groups.iteritems():
        yield pick_best([os.path.join(directory, f) for f in imap(output_filename, repeat(fragment), repeat(cross_fold_index), seeds)])


def make_logo_for_glam2_output(filename):
    "Writes a logo to a filename with .png extension."
    output = GLAM2Output.parse(open(filename))
    freqs, gaps = output.freqs_and_gaps()
    logo = L.pssm_as_image(freqs, size=None, transparencies=gaps)
    logo_filename = '%s.png' % os.path.splitext(filename)[0]
    logo.save(logo_filename)


if '__main__' == __name__:
    for best in best_files('GLAM2'):
        logging.info(best)
        make_logo_for_glam2_output(best)

    from count_gapped_pwm_sites import *
    output = GLAM2Output.parse(open(best))
    freqs, gaps = output.freqs_and_gaps()
    model = build_hmm_model(freqs, gaps, .001)
    fragment, cross_fold_index, seed = interpret_output_filename(os.path.basename(best))
    positive_seqs_filename = os.path.join('Data', '%strimRM-test-x%d.fa' % (fragment, cross_fold_index))
    positive_seqs, tally = load_seqs(positive_seqs_filename)
    positive_scores = test_hmm_forward_backward(model, positive_seqs.values())
    negative_seqs_filename = os.path.join('Data', '%strimRM-test-x%d.fa' % ('T99006', cross_fold_index))
    negative_seqs, tally = load_seqs(negative_seqs_filename)
    negative_scores = test_hmm_forward_backward(model, negative_seqs.values())
    #negative_seqs_filename = os.path.join('Data', '%strimRM-test-x%d.fa' % (fragment, cross_fold_index))
    #negative_scores = test_hmm_forward_backward(model, negative_seqs.values())
    roc_points = roc.picked_rocs_from_thresholds(positive_scores, negative_scores)
    P.close('all')
    P.figure()
    roc.plot_roc_points(roc_points, label=make_tag(fragment, cross_fold_index, seed), color='blue', marker='s')
    roc.plot_random_classifier(label='Random')
    roc.label_plot()
    P.legend(loc='lower right')
#    P.savefig('ROC.eps')
#    P.savefig('ROC.png')
    P.show()
