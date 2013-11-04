#!/usr/local/bin/python
#
# Copyright John Reid 2008, 2009
#

"""
Code to run the single gap algorithm
"""

from optparse import OptionParser
import logging, os, sys
from gapped_pssms.single_gap_algorithm import SingleGapAlgorithm, add_algorithm_options
from gapped_pssms.sequence import convert_fasta_sequences
from output import make_output_dir

#
# Set up the logging
#
logging.basicConfig(level=logging.INFO)


#
# Parse the options
#
option_parser = OptionParser()
option_parser.add_option(
  "-f",
  "--fasta",
  dest="fasta",
  help="The fasta file containing the sequences to run."
)
add_algorithm_options(option_parser)
options, args = option_parser.parse_args()

if not options.output_dir:
    raise ValueError('No output directory specified')
if not os.path.exists(options.output_dir):
    os.makedirs(options.output_dir)
log_filename = os.path.join(options.output_dir, '%s.log' % options.tag)
file_handler = logging.FileHandler(log_filename)
file_handler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
logging.getLogger('').addHandler(file_handler)
logging.info('Writing log to %s', log_filename)
logging.info('Command line: %s', ' '.join(sys.argv))
for option in option_parser.option_list:
    if option.dest:
        logging.info('%32s: %s', option.dest, str(getattr(options, option.dest)))

logging.info('Reading sequences from: %s' % options.fasta)
sequences = convert_fasta_sequences(options.fasta)
#sequences = [s[:200] for s in sequences[:10]]


# Run the algorithm
algorithm = SingleGapAlgorithm(options)
algorithm(sequences)
