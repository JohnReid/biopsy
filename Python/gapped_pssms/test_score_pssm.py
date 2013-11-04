#
# Copyright John Reid 2009
#

"""
Code to score GLAM2 models on the test harness.
"""

import logging, sys
from test_harness_2 import TestHarness, build_hmm_model
from glam2 import GLAM2Output
from parse_gapped_format import parse_models, build_hmm_from_semi_parsed
from optparse import OptionParser


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
  "-q",
  dest="quiet",
  action='store_true',
  default=False,
  help="Run quietly."
)
option_parser.add_option(
  "--glam2-format",
  dest="glam2_format",
  default=False,
  action="store_true",
  help="Algorithm output files are in GLAM2 format."
)
options, args = option_parser.parse_args()
if options.quiet:
    logging.getLogger().setLevel(logging.WARNING)
for option in option_parser.option_list:
    if option.dest:
        logging.info('%30s: %30s (%s)', option.dest, str(getattr(options, option.dest)), option.help)


#
# Parse arguments
#
method, pssm_file, fragment, fold = args[:4]
fold = int(fold)
if 5 == len(args):
    bg = args[4]
    dataset = (fragment, fold, bg)
elif 4 == len(args):
    dataset = (fragment, fold)
else:
    print >> sys.stderr, 'USAGE: %s <method> <pssm-file> <fragment> <fold> [<background>]' % sys.argv[0]
    sys.exit(-1)


#
# Set up the test harness
#
harness = TestHarness(options)


#
# Build the model
#
logging.info('Parsing: %s' % pssm_file)
if options.glam2_format:
    output = GLAM2Output.parse(open(pssm_file))
    freqs, gaps = output.freqs_and_gaps()
    model = build_hmm_model(freqs, gaps)
else:
    semi_parsed_models = list(parse_models(open(pssm_file)))
    if len(semi_parsed_models) > 1:
        print >> sys.stderr, 'For the moment we can only handle one model at a time.'
        sys.exit(-1)
    parsed = semi_parsed_models[0]
    logging.info(str(parsed))
    model, traits = build_hmm_from_semi_parsed(parsed)


#
# Run the model
#
harness.run_method_on_dataset(dataset, method, model)
