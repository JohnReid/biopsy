#!/usr/local/bin/python
#
# Copyright John Reid 2008, 2009
#

"""
Code to take a directory of PSSMs and logs and create a latex report.
"""

from gapped_pssms.score_pssms import PssmSet
import glob, logging, re, os, sys, shutil
from optparse import OptionParser
from itertools import imap

logging.basicConfig(level=logging.INFO)

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
  help="# of PSSMs to report at most"
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
option_parser.add_option(
  "-l",
  "--latex-filename",
  dest="latex_filename",
  default='pssm-report.tex',
  help="Filename to write latex output to."
)

#sys.argv = "latex_report.py -p /home/reid/Analysis/GappedPssms/fragments/revised-fragments/output/2009-02-07--19-39-58".split()
options, args = option_parser.parse_args()
for option in option_parser.option_list:
    if option.dest:
        logging.info('%32s: %s', option.dest, str(getattr(options, option.dest)))


if not options.output_dir:
    options.output_dir = os.path.join(options.pssm_dir, 'report-%d' % options.num_pssms)
if not os.path.exists(options.output_dir):
    os.makedirs(options.output_dir)
logging.getLogger('').addHandler(logging.FileHandler(os.path.join(options.output_dir, 'latex-report.log')))


full_latex_filename = os.path.join(options.output_dir, options.latex_filename)
logging.info('Writing latex report to %s' % full_latex_filename)
latex = open(full_latex_filename, 'w')
print >> latex, """
\\newcommand{\\pssm}[4]{
 \\noindent
 \\includegraphics[width=\\textwidth]{#1} \\*
 #1: \\# seqs with sites: $\\frac{#2}{#3}$, \\# sites: #4
 \\vspace{5pt}
}
"""

tag_re = re.compile('(.*).log')
def get_tag_name(match):
    return match.group(1)
tags = map(
    get_tag_name,
    filter(
        None,
        imap(
            tag_re.search,
            imap(
                os.path.basename,
                glob.glob(os.path.join(options.pssm_dir, '*.log'))
            )
        )
    )
)
for tag in tags:
    print >> latex, '\\section{%s}' % tag
    logging.info(tag)
    pssm_set = PssmSet(os.path.join(options.pssm_dir, tag))
    sorted_pssms = pssm_set.sorted_by_score()
    for i, pssm in enumerate(sorted_pssms[:options.num_pssms]):
        print >> latex, '\\pssm{%s-%03d}{%d}{%d}{%d}\n' % (tag, pssm.pssm_idx, pssm.num_seqs_with_site, pssm.num_seqs, pssm.num_sites)
        shutil.copy(pssm.eps_file, options.output_dir)
        shutil.copy(pssm.png_file, options.output_dir)
latex.close()
