#
# Copyright John Reid 2008, 2009
#

"""
Code to take a directory of PSSMs and logs and rescore the PSSMs, copying the best PSSMs to another directory
for use in the test harness.
"""

from gapped_pssms.score_pssms import PssmSet
import glob, logging, re, os, sys
from optparse import OptionParser
from itertools import imap

#
# Parse the options
#
option_parser = OptionParser()
option_parser.add_option(
  "-n",
  "--num_pssms",
  dest="num_pssms",
  default=5,
  type='int',
  help="# of PSSMs to rescore at most"
)
option_parser.add_option(
  "-p",
  "--pssm-dir",
  dest="pssm_dir",
  default=None,
  help="Directory where the PSSMs are stored"
)
option_parser.add_option(
  "-o",
  "--output-dir",
  dest="output_dir",
  default=None,
  help="Directory where output is written"
)


#
# Parse options and set logs up
#
options, args = option_parser.parse_args()
logging.basicConfig(level=logging.INFO)
if not options.output_dir:
    options.output_dir = os.path.join(options.pssm_dir, 'rescored-%d' % options.num_pssms)
if not os.path.exists(options.output_dir):
    os.makedirs(options.output_dir)
log_filename = os.path.join(options.output_dir, 'rescore-pssms.log')
file_handler = logging.FileHandler(log_filename)
file_handler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
logging.getLogger('').addHandler(file_handler)
logging.info('Command line: %s', ' '.join(sys.argv))
for option in option_parser.option_list:
    if option.dest:
        logging.info('%32s: %s', option.dest, str(getattr(options, option.dest)))


tag_re = re.compile('([TH].*).log')
tags = map(lambda m: m.group(1), filter(None, imap(tag_re.search, glob.glob(os.path.join(options.pssm_dir, '[TH]*.log')))))
for tag in tags:
    logging.info('Looking for PSSMs in %s tagged with %s', options.pssm_dir, tag)
    pssm_set = PssmSet(os.path.join(options.pssm_dir, tag))
    #for pssm in pssm_set.pssms.values():
    #  pssm.write_image()
    sorted_pssms = pssm_set.sorted_by_score()
    output_file = open(os.path.join(options.output_dir, '%s.pssm' % tag), 'w')
    for pssm in sorted_pssms[:options.num_pssms]:
        output_file.write(open(pssm.pssm_file).read())
    output_file.close()
