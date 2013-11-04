#
# Copyright John Reid 2008, 2009
#


"""
A test harness to evaluate gapped pssm models.
"""


import os, os.path, sys, logging, numpy, glob, random, gp_site_config, cookbook, cPickle
from optparse import OptionParser
from time import strftime
from infpy.roc import RocCalculator, update_roc, get_new_roc_parameter, plot_roc_points, label_plot, plot_random_classifier
from hmm.pssm import seq_to_numpy
from gapped_pssms.parse_gapped_format import parse_models, build_hmm_from_semi_parsed
from itertools import imap, islice, izip


class DictOfRocs(dict):
    """A dictionary of ROCs"""
    def __missing__(self, k):
        self[k] = RocCalculator()
        return self[k]



def ensure_dir_exists(dir):
    if os.access(dir, os.X_OK):
        return
    else:
        # make sure parent exists
        parent_dir = os.path.dirname(dir)
        if parent_dir:
            ensure_dir_exists(parent_dir)
            os.mkdir(dir)





def evaluate_model(model, sequence):
    """
    Evaluates the model against the sequence.

    @return: True if there is at least one hit in the sequence
    """
    hmm, traits = model
    LL, states = hmm.viterbi(seq_to_numpy(sequence))
    # we have a hit if we find at least K/2 states in the state sequence that are not in the
    # background
    return sum(state not in traits.background_states for state in states) > traits.K / 2



def evaluate_models(models, sequence):
    """
    Evaluates the models against the sequence.

    @return: The index of the first model to have a hit in the sequence or None
    """
    for i, model in enumerate(models):
        if evaluate_model(model, sequence):
            return i
    return None



def calculate_model_counts_on_sequences(models, sequence_set):
    """
    Evaluate the models on the sequences. Make a count of how many sequences each model was the first to find a hit in.
    """
    result = cookbook.DictOfInts()
    for s in sequence_set:
        first_pssm = evaluate_models(models, s)
        if None != first_pssm:
            result[first_pssm] += 1
    return result






def score_sequence_set(sequence_set, pssms, p_binding_site):
    models = [
        build_hmm_from_semi_parsed(
            parsed,
            p_binding_site=p_binding_site
        )
        for parsed
        in pssms
    ]
    return calculate_model_counts_on_sequences(models, sequence_set)



def score_sequences(sequences, pssms, bg_types, p_binding_sites):
    """
    Take a bunch of PSSMs and score sequence sets with them.
    """
    scores = cookbook.DictOfLists()
    for (dataset, cross_fold_index), pssms_for_fold in pssms.iteritems():
        positive_key = positive_sequences_key(dataset, cross_fold_index)
        positive_sequence_set = sequences[positive_key]
        logging.info(
          'Analysing %d positive sequences (%d bases) for %s-%d',
          len(positive_sequence_set),
          sum(len(s) for s in positive_sequence_set),
          dataset,
          cross_fold_index
        )
        for p_binding_site in p_binding_sites:
            scores[positive_key].append(score_sequence_set(positive_sequence_set, pssms_for_fold, p_binding_site))

        for bg_type in bg_types:
            negative_key = negative_sequences_key(dataset, cross_fold_index, bg_type)
            negative_sequence_set = sequences[negative_key]
            logging.info(
              'Analysing %d negative sequences (%d bases) for %s-%d-%s',
              len(negative_sequence_set),
              sum(len(s) for s in negative_sequence_set),
              dataset,
              cross_fold_index,
              bg_type
            )
            for p_binding_site in p_binding_sites:
                scores[negative_key].append(score_sequence_set(negative_sequence_set, pssms_for_fold, p_binding_site))

    return scores



def score_to_count(score, num_pssms):
    """
    A score is a dict. The keys index PSSMs, the values are the counts of when that PSSM was the first one to have a hit in a given sequence.
    We would like to have counts of when that PSSM or one with a lower index had a hit. This function does that conversion.
    """
    result = [score[i] for i in xrange(num_pssms)]
    for i in xrange(1, num_pssms):
        result[i] += result[i-1]
    return result


def scores_to_counts(scores, num_pssms):
    return dict(
        (key, [score_to_count(score, num_pssms) for score in value])
        for key, value in scores.iteritems()
    )



def sequences_from_fasta(fasta):
    """Yields sequences from fasta file."""
    import corebio.seq_io.fasta_io
    return imap(
      lambda s: s.strip('nN'),
      imap(
        str,
        corebio.seq_io.fasta_io.iterseq(
          open(fasta, 'r'),
          corebio.seq.dna_alphabet
        )
      )
    )

def positive_sequence_filename(dataset, cross_fold_index):
    "@return: The filename for the positive sequences for the given fragment and cross fold index."
    return os.path.join(
        options.positive_dir,
        options.positive_sequence_pattern % (dataset, cross_fold_index+1)
    )

def bg_sequence_filename(bg_type):
    "@return: The filename for the background sequences for the given type."
    return os.path.join(
        options.bg_dir,
        options.bg_sequence_pattern % bg_type
    )

def sequence_generator_for_bg_type(bg_type):
    "Yields sequences for this background type."
    filename = bg_sequence_filename(bg_type)
    logging.info('Loading background sequences for %s from %s', bg_type, filename)
    while True:
        for s in sequences_from_fasta(filename):
            yield s

def positive_sequences_key(dataset, cross_fold_index):
    "@return: A key for indexing a set of positive sequences for this fragment and cross fold index."
    return "positive-%s-%d" % (dataset, cross_fold_index)

def negative_sequences_key(dataset, cross_fold_index, bg_type):
    "@return: A key for indexing a set of negative sequences for this fragment and cross fold index and given background type."
    return "negative-%s-%d-%s" % (dataset, cross_fold_index, bg_type)

def assemble_sequences(datasets, num_folds, bg_types):
    """
    Assemble a dict mapping keys to sets of sequences for positive and negative examples.
    """

    # the result
    sequences = dict()

    # get sequences for all background types
    bg_sequence_generators = dict((bg_type, sequence_generator_for_bg_type(bg_type)) for bg_type in bg_types)

    # for each dataset
    for dataset in datasets:

        # for each cross fold
        for cross_fold_index in xrange(num_folds):

            # get the positive sequences
            positive_sequences = list(sequences_from_fasta(positive_sequence_filename(dataset, cross_fold_index)))
            sequences[positive_sequences_key(dataset, cross_fold_index)] = positive_sequences

            # for each negative type
            for bg_type in bg_types:

                def shorten_negative(n, p):
                    if len(n) < len(p):
                        logging.warning('Negative sequence too short to match positive sequence of length: %d; dataset: %s-%d', len(p), dataset, cross_fold_index)
                        # raise RuntimeError('Negative sequence too short to match positive sequence of length: %d' % len(p))
                    return n[:len(p)]

                # match the length of each sequence to the lengths of the positive sequences
                negative_sequences = [shorten_negative(n, p) for n, p in izip(bg_sequence_generators[bg_type], positive_sequences)]
                if len(negative_sequences) < len(positive_sequences):
                    raise RuntimeError('Not enough negative sequences')

                # store the sequences in our dict
                sequences[negative_sequences_key(dataset, cross_fold_index, bg_type)] = negative_sequences

    return sequences



def look_for_pssms_in_dir(directory):
    """Find the models that are in the specified directory (named <dataset>-<cross fold index>.pssm)."""
    # for each file called *.pssm in the model directory
    logging.info('Looking for PSSMs in: %s' % directory)
    parsed_models = dict()
    for file in glob.glob(os.path.join(directory, '*-*.pssm')):
        try:
            #logging.info('Found models in %s' % file)
            base, ext = os.path.splitext(os.path.basename(file))
            dataset, cross_fold_idx = base.split('-')
            cross_fold_idx = int(cross_fold_idx)-1
            logging.info(
              'Read models from: %s for dataset: %s and cross-fold validation set: %d',
              file,
              dataset,
              cross_fold_idx
            )
            parsed_models[(dataset, cross_fold_idx)] = list(
              parse_models(
                open(os.path.join(file))
              )
            )
        except:
            print sys.exc_info()
            logging.warning('Could not parse: %s' % file)
    return parsed_models




def calculate_p_binding_site_parameters(num_roc_points, min_p_binding_site, max_p_binding_site):
    """Calculate a range of values for the p_binding_site parameter."""
    log_min_p_binding_site = numpy.log10(min_p_binding_site)
    log_max_p_binding_site = numpy.log10(max_p_binding_site)
    return numpy.power(
      10,
      numpy.arange(
        log_min_p_binding_site,
        log_max_p_binding_site,
        (log_max_p_binding_site-log_min_p_binding_site)/num_roc_points
      )
    )[:num_roc_points]




def calculate_rocs(counts, bg_types, datasets, num_points, num_folds, num_pssms):
    """
    Calculate ROC curves.
    """
    rocs = dict()
    for bg_type in bg_types:
        rocs[bg_type] = dict()
        for dataset in datasets:
            rocs[bg_type][dataset] = dict()
            for p in xrange(num_pssms):
                rocs[bg_type][dataset][p] = roc_points = [RocCalculator() for i in xrange(num_points)]
                for cross_fold_index in xrange(num_folds):
                    # add the positive test examples
                    positive_key = positive_sequences_key(dataset, cross_fold_index)
                    N = len(sequences[positive_key])
                    assert len(counts[positive_key]) == num_points
                    for rp, c in izip(roc_points, counts[positive_key]):
                        rp.tp += c[p]
                        rp.fn += N-c[p]
                    # add the positive test examples
                    negative_key = negative_sequences_key(dataset, cross_fold_index, bg_type)
                    assert len(counts[negative_key]) == num_points
                    for rp, c in izip(roc_points, counts[negative_key]):
                        rp.fp += c[p]
                        rp.tn += N-c[p]
    return rocs


def plot_rocs(title, rocs, num_pssms):
    """
    Plot roc graphs.
    """
    import pylab as P
    for bg_type, rocs_for_bg in rocs.iteritems():
        for dataset, rocs_for_dataset in rocs_for_bg.iteritems():
            tag = '%s-%s-%s' % (title, bg_type, dataset)
            fig = P.figure()
            for p in num_pssms:
                plot_roc_points(rocs_for_dataset[p][::-1], label='%d PSSMs' % (p+1))
            P.legend(loc='lower right')
            P.title(tag)
            plot_random_classifier()
            label_plot()
            P.savefig(os.path.join(options.output_dir, 'roc-%s.png' % tag))
            P.close(fig)




def write_roc_statistics(rocs, p_binding_sites, num_pssms):
    """
    Write the statistics for the ROCs.
    """
    filename = os.path.join(options.output_dir, options.roc_statistics_file)
    logging.info('Writing ROC statistics to %s', filename)
    f = open(filename, 'w')
    f.write('# bg-type; dataset; num pssms;p(binding);TP;TN;FP;FN\n')
    for bg_type, rocs_for_bg in rocs.iteritems():
        for dataset, rocs_for_dataset in rocs_for_bg.iteritems():
            tag = '%s-%s' % (bg_type, dataset)
            for p in num_pssms:
                assert len(rocs_for_dataset[p]) == len(rocs_for_dataset[p])
                for roc, p_binding_site in zip(rocs_for_dataset[p], p_binding_sites):
                    f.write('%s;%s;%d;%f;%d;%d;%d;%d\n' % (bg_type, dataset, p, p_binding_site, roc.tp, roc.tn, roc.fp, roc.fn))
    f.close()











logging.basicConfig(level=logging.INFO)
logging.getLogger('').addHandler(logging.FileHandler('test_harness.log'))

#
# Parse the options
#
usage = """usage: %prog [options] <bg sequence names>

Put all of your PSSMs in a directory in files with a .pssm extension. The PSSMs for the 3rd fold
for dataset T99001 would be called T99001-3.pssm for example. The directory is specified by the -m
option.

Decide which background models you would like to run ROCs for, e.g.
vm-r1-back-human-hsNCBI36-v2 vm-r3-TSS-human-hsNCBI36-v2
and add these to the end of the command line arguments.

By default the output will be placed in a sub-directory of the directory the models are in.
"""
option_parser = OptionParser(usage=usage)
option_parser.add_option(
  "-m",
  "--model-dir",
  dest="model_dir",
  help="A directory containing the gapped pssm model(s)"
)
option_parser.add_option(
  "-o",
  "--output-dir",
  dest="output_dir",
  default=None,
  help="Directory where output is written"
)
option_parser.add_option(
  "--positive-filename-pattern",
  dest="positive_sequence_pattern",
  default=gp_site_config.default_positive_sequence_test_pattern,
  help="A printf-style format string that when given the dataset and the cross-fold index generates the filename of its sequences"
)
option_parser.add_option(
  "--positive-directory",
  dest="positive_dir",
  default=gp_site_config.default_positive_sequence_dir,
  help="The directory containing the positive sequences"
)
option_parser.add_option(
  "--bg-filename-pattern",
  dest="bg_sequence_pattern",
  default='%s.fa',
  help="A printf-style format string that when given one of the background types generates the filename of its sequences"
)
option_parser.add_option(
  "--bg-directory",
  dest="bg_dir",
  default=gp_site_config.default_background_sequence_dir,
  help="The directory containing the negative sequences"
)
#option_parser.add_option(
#  "-r",
#  "--roc-output-file",
#  default='test-harness-roc',
#  dest="roc_output_file",
#  help="The filename to write the ROC curve to"
#)
option_parser.add_option(
  "-t",
  "--title",
  dest="title",
  help="The title of test harness run."
)
option_parser.add_option(
  "-r",
  "--roc-statistics-file",
  default='roc-statistics.txt',
  dest="roc_statistics_file",
  help="Filename to write ROC statistics to"
)
option_parser.add_option(
  "--num-points",
  dest="num_points",
  default=16,
  type='int',
  help="The number of points to generate on the ROC curve"
)
option_parser.add_option(
  "--max-p-binding-site",
  dest="max_p_binding_site",
  default=0.1,
  type=float,
  help="The largest p(binding site) to use in Viterbi algorithm, will generate right-most point on ROC."
)
option_parser.add_option(
  "--min-p-binding-site",
  dest="min_p_binding_site",
  default=0.0001,
  type=float,
  help="The smallest p(binding site) to use in Viterbi algorithm, will generate left-most point on ROC."
)
#option_parser.add_option(
#  "--random-seed",
#  dest="random_seed",
#  default=None,
#  type=int,
#  help="Seeds the random shuffle of the background sequences."
#)
option_parser.add_option(
  "--num-folds",
  dest="num_folds",
  default=5,
  type=int,
  help="The number of cross-validation folds."
)
#sys.argv='test_harness.py -m /home/john/Analysis/GappedPssms/fragments/x-validation-test vm-r1-back-human-hsNCBI36-v2 vm-r3-TSS-human-hsNCBI36-v2'.split()
option_parser.add_option(
  "--num-pssms",
  dest="num_pssms",
  default=5,
  type=int,
  help="The number of PSSMs to use"
)
#sys.argv='test_harness.py -m /home/john/Analysis/GappedPssms/fragments/x-validation-test vm-r1-back-human-hsNCBI36-v2 vm-r3-TSS-human-hsNCBI36-v2'.split()
options, args = option_parser.parse_args()

#
# Sort out the output directory
#
if not options.output_dir:
    options.output_dir = os.path.join(options.model_dir, 'test-harness')
if not os.path.exists(options.output_dir):
    os.makedirs(options.output_dir)
logging.getLogger('').addHandler(logging.FileHandler(os.path.join(options.output_dir, 'test-harness.log')))


#
# Log what the program options are
#
logging.info('Command line: %s', ' '.join(sys.argv))
for option in option_parser.option_list:
    if option.dest:
        logging.info('%32s: %s', option.dest, str(getattr(options, option.dest)))


#if options.random_seed:
#    logging.info('Seeding randomness with %d', options.random_seed)
#    random.seed(options.random_seed)
#else:
#    logging.info('Using random seed')
#    random.seed()

#
# Find the PSSMs in their directory
#
pssms = look_for_pssms_in_dir(options.model_dir)


#
# Work out what datasets we have PSSMs for and what the background types are
#
datasets = set(dataset for (dataset, cross_fold_index) in pssms.keys())
bg_types = args
if not len(bg_types):
    raise RuntimeError('No background types specified on command line.')


#
# Calculate which p_binding_site parameters to use.
#
p_binding_sites = calculate_p_binding_site_parameters(options.num_points, options.max_p_binding_site, options.min_p_binding_site)
logging.info('Using the following p(binding) values: %s', str(p_binding_sites))


#
# Assemble positive and negative sequence sets
#
sequences = assemble_sequences(datasets, options.num_folds, bg_types)



#
# Score the sequences using the PSSMs
#
scores = score_sequences(sequences, pssms, bg_types, p_binding_sites)
scores_pickle_file = os.path.join(options.output_dir, 'scores.pickle')
logging.info('Storing scores in %s' % scores_pickle_file)
cPickle.dump(scores, open(scores_pickle_file, 'w'))


#
# Convert the scores to counts
#
counts = scores_to_counts(scores, options.num_pssms)



#
# Calculate the ROCs
#
rocs = calculate_rocs(counts, bg_types, datasets, len(p_binding_sites), options.num_folds, options.num_pssms)


#
# Plot the ROCs
#
plot_rocs(options.title, rocs, range(options.num_pssms))



#
# Write the ROC statistics
#
write_roc_statistics(rocs, p_binding_sites, range(options.num_pssms))




logging.info('Test harness finished.')
