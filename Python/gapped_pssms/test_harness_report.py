#
# Copyright John Reid 2009
#


"""
Code to make a report from the test harness.
"""


import os, logging, sys, glam2, test_harness_methods
from optparse import OptionParser
from test_harness_2 import TestHarness

#
# Initialise the logging
#
format='%(asctime)-15s %(message)s'
logging.basicConfig(level=logging.INFO, format=format)
logging.info('Command line: %s', ' '.join(sys.argv))
logging.info('Current working directory: %s', os.getcwd())


#
# Parse the options
#
option_parser = OptionParser()
TestHarness.add_options(option_parser)
option_parser.add_option(
  "--roc-dir",
  dest="roc_dir",
  default='.',
  help="Directory ROC files are stored in."
)
options, args = option_parser.parse_args()
for option in option_parser.option_list:
    if option.dest:
        logging.info('%30s: %30s (%s)', option.dest, str(getattr(options, option.dest)), option.help)
dirs = map(test_harness_methods.load_dir, args)
keys = set()
for dir in dirs:
    keys.update(dir.files.keys())
fragments = set(f for f, i in keys)

tex = open('methods.tex', 'w')
for fragment in fragments:
    print >> tex, """
%s\\\\

\\includegraphics{%s}
""" % (
    fragment,
    os.path.join(options.roc_dir, 'ROC-%s-shuffle' % fragment),
)
tex.close()
