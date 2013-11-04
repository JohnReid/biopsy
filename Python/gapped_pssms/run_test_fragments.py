#
# Copyright John Reid 2008, 2009
#

"""
Code to run the single gap algorithm on test fragments.
"""


from optparse import OptionParser
import logging, os, sys
from gapped_pssms.single_gap_algorithm import SingleGapAlgorithm, add_algorithm_options
from gapped_pssms.sequence import convert_fasta_sequences
from gapped_pssms.data import fasta_file_for_fragment, test_set_fragments

#
# Set up the logging
#
logging.basicConfig(level=logging.DEBUG)


#
# Parse the options
#
option_parser = OptionParser()
add_algorithm_options(option_parser)
logging.info('Command line: %s', ' '.join(sys.argv))
options, args = option_parser.parse_args()

log_filename = os.path.join(options.output_dir, '%s.log' % options.tag)
logging.getLogger('').addHandler(logging.FileHandler(log_filename))
logging.info('Writing log to %s', log_filename)
for option in option_parser.option_list:
    if option.dest:
        logging.info('%s: %s (%s)', option.dest, str(getattr(options, option.dest)), option.help)

#inputs = [('K10-g0.50-N200-L200-seed4-1', fasta_file_for_synthetic_data('K10-g0.50-N200-L200-seed4-1'))]
inputs = [(fragment, fasta_file_for_fragment(fragment)) for fragment in test_set_fragments]

# for each input sequence
for seq_tag, fasta_file in inputs:
    # add a file handler to log for this test set
    file_handler = logging.FileHandler(os.path.join(options.output_dir, '%s.log' % seq_tag))
    logging.getLogger('').addHandler(file_handler)
    try:
        sequences = convert_fasta_sequences(fasta_file)
        #sequences = [s[:200] for s in sequences[:10]]

        # set up the options for this test set
        options.tag = seq_tag
        options.bg_model_filename = "%s-bg-model.pickle" % seq_tag

        # Run the algorithm
        algorithm = SingleGapAlgorithm(options)
        algorithm(sequences)

    finally:
        logging.getLogger('').removeHandler(file_handler)
        file_handler.close()
