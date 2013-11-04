#
# Copyright John Reid 2009
#

"""
Code to implement a test harness.
"""


import corebio.seq_io.fasta_io, os, logging, hmm, cPickle, numpy as N, sys, numpy.random as R
from hmm.pssm import seq_to_numpy, numpy_to_seq
from optparse import OptionParser
from cookbook.dicts import DictOf
from itertools import imap, cycle

_logger=logging.getLogger(__name__)


def fragment_name(fragment):
    if 'T00671' == fragment:
        return 'p53'
    if 'T00759' == fragment:
        return 'Sp1'
    if 'T99002' == fragment:
        return 'GABP'
    if 'T99003' == fragment:
        return 'NRSF'
    if 'T99005' == fragment:
        return 'Stat5a'
    if 'T99006' == fragment:
        return 'Stat5b'
    return fragment


def load_sequences(fasta):
    "Load sequences."
    _logger.info('Loading sequences: %s', fasta)
    sequences = [
        seq.remove('Nn')
        for seq
        in corebio.seq_io.fasta_io.iterseq(open(fasta, 'r'), corebio.seq.reduced_nucleic_alphabet)
    ]
    return sequences


def is_unknown(x):
    "@return: True iff x == 4."
    return 4 == x


def partition(iterable, keyfunc=None):
    "Partition an iterable into groups for which keyfunc returns the same value. Yields (key, begin, end) tuples."
    if keyfunc is None:
        keyfunc = lambda x: x
    iterable = iter(iterable)
    lastkey = object()
    begin = None
    for i, x in enumerate(iterable):
        currkey = keyfunc(x)
        if currkey != lastkey:
            if None != begin:
                yield lastkey, begin, i
            begin = i
        lastkey = currkey
    if None != begin:
        yield lastkey, begin, i


def shuffle_sequence(seq):
    "@return: A shuffled version of the sequence leaving unknowns in place."
    ords = N.array(seq.ords())
    for key, begin, end in partition(ords, keyfunc=is_unknown):
        if not key:
            R.shuffle(ords[begin:end])
    shuffled = seq.alphabet.chrs(ords)
    shuffled.name = '%s (shuffled)' % seq.name
    shuffled.description = '%s (shuffled)' % seq.description
    return shuffled


def shorten_seq(to_shorten, sequence_to_match_length):
    "Shorten a sequence to match the length of the other."
    if len(to_shorten) < len(sequence_to_match_length):
        raise RuntimeError('Not enough bases in sequence: %d < %d' % (len(to_shorten), len(sequence_to_match_length)))
    ords = N.array(to_shorten.ords()[:len(sequence_to_match_length)])
    ords[N.where(4==N.array(sequence_to_match_length.ords()))] = 4 # set unknown bases same as original sequence
    result = sequence_to_match_length.alphabet.chrs(ords)
    result.name = '%s (matched)' % sequence_to_match_length.name
    result.description = '%s (matched)' % sequence_to_match_length.description
    return result



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
        gamma0[N.where(sequence.as_numpy()==4)] = 1. # make sure that where we have Ns aren't counted
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


class TestHarness(object):
    "A test harness."


    def __init__(self, options):
        "Construct."

        self.options = options
        "Options for the test harness."

        # put defaults in if not specified
        if not len(options.fragments):
            self.options.fragments = default_fragments
        if not len(options.backgrounds):
            self.options.backgrounds = default_backgrounds

        self.lazy_sequences = DictOf(load_sequences, take_key_as_arg=True)
        "Reads sequences lazily."

        self.lazy_sequences_for_hmm = DictOf(self.get_sequence_for_hmm, take_key_as_arg=True)
        "Converts sequences to HMM format lazily."


    def get_sequence_for_hmm(self, fasta):
        "Get sequences converted to how HMM module likes them."
        _logger.info('Converting sequences: %s', fasta)
        sequences = self.lazy_sequences[fasta]
        numpy_seqs = map(hmm.preprocess_sequence, imap(N.array, imap(corebio.seq.Seq.ords, sequences)))
        return numpy_seqs


    def tag(self, dataset):
        "@return: A tag for the dataset."
        if 2 == len(dataset):
            return '%s-x%d' % dataset # positive data set
        elif 3 == len(dataset):
            return '%s-x%d-neg-%s' % dataset # negative data set
        else:
            raise RuntimeError('Unexpected length of data set arguments')


    def fasta_filename(self, dataset):
        "@return: The fasta filename for the data set."
        return os.path.join(self.options.data_dir, '%s.fa' % self.tag(dataset))


    def sequences(self, dataset):
        "@return: The sequences for the dataset."
        return self.lazy_sequences[self.fasta_filename(dataset)]


    def hmm_sequences(self, dataset):
        "@return: The HMM sequences for the dataset."
        return self.lazy_sequences_for_hmm[self.fasta_filename(dataset)]


    def results_filename(self, dataset, method):
        "@return: The results filename."
        return os.path.join(self.options.results_dir, '%s-%s.results' % (self.tag(dataset), method))


    def results(self, dataset, method):
        "@return: Load results from disk."
        return cPickle.load(open(self.results_filename(dataset, method), 'rb'))


    def run_method_on_dataset(self, dataset, method, model):
        "Run the method on the dataset."
        logging.info('Running %s on %s', method, dataset)
        seqs = self.hmm_sequences(dataset)
        scores = test_hmm_forward_backward(model, seqs)
        filename = self.results_filename(dataset, method)
        logging.info('First results are: %s', scores[:4])
        logging.info('Writing results to %s', filename)
        cPickle.dump(scores, open(filename, 'wb'), protocol=cPickle.HIGHEST_PROTOCOL)


    def all_datasets(self):
        "@return: Yield all the datasets."
        for fragment in self.options.fragments:
            for fold in self.folds():
                yield (fragment, fold)
                for bg in self.options.backgrounds:
                    yield  (fragment, fold, bg)


    def folds(self):
        "@return: The folds."
        return range(1, self.options.num_folds+1)


    def build_data_dir(self, original_dir, backgrounds):
        """
        Create and populate a directory of fasta files for the test harness.

        Link positive test fasta files to originals.

        Read background fasta files in and create negative test sets of correct lengths.
        """
        # make sure data directory exists
        if not os.path.exists(self.options.data_dir):
            os.makedirs(self.options.data_dir)

        # create a dictionary of generators, one for each background
        background_seqs = dict(
            (bg, cycle(corebio.seq_io.fasta_io.iterseq(open(fasta), corebio.seq.reduced_nucleic_alphabet)))
            for bg, fasta
            in backgrounds.iteritems()
        )

        def get_output_file(dataset):
            "Get output filename for dataset, removing existing file if it exists."
            filename = self.fasta_filename(dataset)
            if os.path.exists(filename):
                os.remove(filename)
            return filename

        # for each fragment and fold
        for fragment in self.options.fragments:
            for fold in self.folds():

                # link positive sequences to original sequences
                dataset = (fragment, fold)
                original_fasta = os.path.join(original_dir, '%strimRM-test-x%d.fa' % dataset)
                os.symlink(original_fasta, get_output_file(dataset))
                positive_seqs = self.sequences(dataset)


                # get a set of negative sequences for each background and write to new file
                for bg, bg_seqs in background_seqs.iteritems():
                    dataset = (fragment, fold, bg)
                    bg_seqs = corebio.seq.SeqList(
                        imap(shorten_seq, bg_seqs, positive_seqs),
                        name=None,
                        description=None,
                        alphabet=corebio.seq.reduced_nucleic_alphabet,
                    )
                    if len(bg_seqs) != len(positive_seqs):
                        raise RuntimeError('Not enough background negative sequences %s: %d != %d' % (bg, len(bg_seqs), len(positive_seqs)))
                    corebio.seq_io.fasta_io.write(open(get_output_file(dataset), 'w'), bg_seqs)

                # create a negative dataset of shuffled versions of the positive dataset
                dataset = (fragment, fold, 'shuffle')
                f = open(get_output_file(dataset), 'w')
                for seq in positive_seqs:
                    corebio.seq_io.fasta_io.writeseq(f, shuffle_sequence(seq))
                f.close()


    @staticmethod
    def add_options(option_parser):
        "Add test harness options to the parser."
        option_parser.add_option(
            '-d',
            '--data-dir',
            dest="data_dir",
            default='.',
            help="Where the sequences are stored."
        )
        option_parser.add_option(
            '-r',
            '--results-dir',
            dest="results_dir",
            default='.',
            help="Where the test harness results are stored."
        )
        option_parser.add_option(
            '-b',
            dest="backgrounds",
            default=[],
            action='append',
            help="The background data sets. If none specified defaults are used: %s" % default_backgrounds
        )
        option_parser.add_option(
            '-f',
            dest="fragments",
            default=[],
            action='append',
            help="The fragments. If none specified defaults are used: %s" % default_fragments
        )
        option_parser.add_option(
            '-n',
            '--num-folds',
            dest="num_folds",
            type='int',
            default=5,
            help="How many folds in the cross validation."
        )


default_fragments = [
    'T00671',
    'T00759',
    'T99002',
    'T99003',
    'T99005',
    'T99006',
]

default_backgrounds = [
    'r1-back',
    'r3-TSS',
    'shuffle',
]

def choose_existing_dir(candidate_dirs):
    "Given a list of candidate directories, chooses the first one that exists or returns None."
    for dir in candidate_dirs:
        if os.path.exists(dir):
            return dir
    return None


if '__main__' == __name__:
    #
    # Initialise the logging
    #
    logging.basicConfig(level=logging.INFO)
    log_filename = 'create-test-harness-data.log'
    file_handler = logging.FileHandler(log_filename)
    file_handler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
    logging.getLogger('').addHandler(file_handler)
    logging.info('Writing log to %s', log_filename)
    logging.info('Command line: %s', ' '.join(sys.argv))

    data_dir = choose_existing_dir(
        [
            '/home/reid/Data/GappedPssms',
            '/home/john/Data/GappedPssms',
        ]
    )

    option_parser = OptionParser()
    option_parser.add_option(
        '--background-directory',
        dest="bg_dir",
        default=os.path.join(data_dir, 'neg-sequences'),
        help="The directory with the background (negative sequences)."
    )
    option_parser.add_option(
        '--original-directory',
        dest="orig_dir",
        default=os.path.join(data_dir, 'apr-2009'),
        help="The directory with original sequences."
    )
    TestHarness.add_options(option_parser)
    options, args = option_parser.parse_args()
    for option in option_parser.option_list:
        if option.dest:
            logging.info('%32s: %s', option.dest, str(getattr(options, option.dest)))

    harness = TestHarness(options)
    backgrounds = {
        'r1-back' : os.path.join(options.bg_dir, 'vm-r1-back-human-hsNCBI36-v2.fa'),
        'r3-TSS'  : os.path.join(options.bg_dir, 'vm-r3-TSS-human-hsNCBI36-v2.fa'),
    }
    harness.build_data_dir(options.orig_dir, backgrounds)
