#
# Copyright John Reid 2008, 2009
#

"""
Applies single gap motif finding algorithm to cross-validation training data sets.
"""

from optparse import OptionParser
import logging, os
from gapped_pssms.single_gap_algorithm import SingleGapAlgorithm, add_algorithm_options
from gapped_pssms.sequence import convert_fasta_sequences
from gapped_pssms.data import fasta_files_for_fragment_cross_fold, test_set_fragments, cross_folds

#
# Set up the logging
#
logging.basicConfig(level=logging.INFO)


#
# Parse the options
#
option_parser = OptionParser()
add_algorithm_options(option_parser)
options, args = option_parser.parse_args()

if not options.output_dir:
    options.output_dir = make_output_dir(os.path.join(results_dir, 'single-gap')),
if not os.path.exists(options.output_dir):
    os.makedirs(options.output_dir)


log_filename = os.path.join(options.output_dir, '%s.log' % options.tag)
file_handler = logging.FileHandler(log_filename)
file_handler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
logging.getLogger('').addHandler(file_handler)
logging.info('Command line: %s', ' '.join(sys.argv))
logging.info('Writing log to %s', log_filename)
for option in option_parser.option_list:
    if option.dest:
        logging.info('%32s: %s', option.dest, str(getattr(options, option.dest)))


#
# For each data set
#
for fragment in test_set_fragments:
    for i in xrange(1, cross_folds+1):
        train_fasta, test_fasta = fasta_files_for_fragment_cross_fold(fragment, i)
        logging.info('Training data set: %s', train_fasta)
        #log.info('Test data set: %s', test_fasta)

        seq_tag = '%s-%d' % (fragment, i)
        file_handler = logging.FileHandler(os.path.join(options.output_dir, '%s.log' % seq_tag))
        logging.getLogger('').addHandler(file_handler)
        try:
            sequences = convert_fasta_sequences(train_fasta)
            #sequences = [s[:200] for s in sequences[:100]]

            # set up the options for this test set
            options.tag = seq_tag
            options.bg_model_filename = "%s-bg-model.pickle" % seq_tag

            # Run the algorithm
            algorithm = SingleGapAlgorithm(options)
            algorithm(sequences)

        finally:
            logging.getLogger('').removeHandler(file_handler)
            file_handler.close()
