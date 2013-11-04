#
# Copyright John Reid 2008
#


"""
Code to parse single gapped PSSMs from a file and write their logos out.
"""

import os, sys, numpy, logging
from optparse import OptionParser
from gapped_pssms.parse_gapped_format import parse_models, build_hmm_from_semi_parsed
import hmm.pssm.logo as logo


def emissions_and_gaps_from_semi_parsed(p):
    emissions = numpy.zeros((p.N, p.M))
    gap_probs = numpy.ones(p.N)
    for i, e in p.emissions:
        emissions[i] = e
    for i, j, t in p.transitions:
        if i == j-1:
            gap_probs[j] = t
    return emissions, gap_probs

#
# Initialise the logging
#
logging.basicConfig(level=logging.INFO)




#
# Parse the options
#
option_parser = OptionParser()
option_parser.add_option(
  "-m",
  "--models-file",
  dest="models_file",
  help="File in which the gapped PSSMs are stored."
)
option_parser.add_option(
  "-l",
  "--logo-files-basename",
  dest="logo_files_basename",
  help="basename of files to write logos to. Extension will be -0.png"
)
option_parser.add_option(
  "-t",
  "--image-type",
  dest="image_type",
  default='png',
  help="type of images to write"
)
options, args = option_parser.parse_args()
for option in option_parser.option_list:
    if option.dest:
        logging.info('%s: %s (%s)', option.dest, str(getattr(options, option.dest)), option.help)


# Load PSSMs
logging.info('Loading PSSMs: %s', options.models_file)
pssms = list(parse_models(open(options.models_file)))

for i, p in enumerate(pssms):
    filename = '%s-%d.%s' % (options.logo_files_basename, i, options.image_type)
    logging.info('Creating image for PSSM: %s', filename)

    emissions, gap_probs = emissions_and_gaps_from_semi_parsed(p)
    logo_image = logo.pssm_as_image(
      emissions,
      transparencies=gap_probs
    )

    logo_image.save(filename)
