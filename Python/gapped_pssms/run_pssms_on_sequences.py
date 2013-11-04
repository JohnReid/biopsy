#
# Copyright John Reid 2008
#


"""
Code to parse single gapped PSSMs from a file and apply them to sequences in a FASTA file.
"""

import logging, os, sys, itertools
from optparse import OptionParser
from gapped_pssms.parse_gapped_format import parse_models, build_hmm_from_semi_parsed
from hmm.pssm import seq_to_numpy, numpy_to_seq
import pylab as P


def sequences_from_fasta(fasta):
    """Yields sequences from fasta file."""
    import corebio.seq_io.fasta_io
    for seq in corebio.seq_io.fasta_io.iterseq(
      open(fasta, 'r'),
      corebio.seq.dna_alphabet
    ):
        yield seq.description.strip(), str(seq).strip('nN')

def initial_bs_states(hmm, traits):
    result = set()
    for s in traits.background_states:
        result.update(set(s for s in numpy.where(hmm.A[s] > 0.)[0]))
    result.difference_update(traits.background_states)
    return result


def analyse_sequence_for_best_site(sequence, hmm, traits):
    LL, alpha, beta, c = hmm.forward_backward(sequence)
    ab = alpha * beta
    p_state = (ab.T / ab.sum(axis=1)).T
    p_not_bg = 1. - p_state[:,traits.background_states].sum(axis=1)
    best_base = p_not_bg.argmax()
    best_score = p_not_bg[best_base]
    start = best_base
    while start >= 0 and p_not_bg[start-1] > .5 * best_score:
        start -= 1
    end = best_base
    while end < len(sequence) and p_not_bg[end] > .5 * best_score:
        end += 1
    return best_score, sequence[start:end]

class SiteFinder(object):
    def __init__(self, hmm, traits, threshold=0.9):
        self.threshold = threshold
        self.hmm = hmm
        self.A = hmm.A
        self.traits = traits

    def __call__(self, sequence):
        self.sequence = sequence
        self._analyse_sequence()
        for site in self._look_for_sites():
            yield site

    def _analyse_sequence(self):
        LL, self.alpha, self.beta, c = self.hmm.forward_backward(self.sequence)
        self.ab = self.alpha * self.beta
        p_state = (self.ab.T / self.ab.sum(axis=1)).T
        self.p_not_bg = 1. - p_state[:,self.traits.background_states].sum(axis=1)
        print self.p_not_bg.max()

    def _look_for_sites(self):
        self.current_base = 0
        while not self._at_end():
            if self._test_for_site_start():
                yield self._calculate_site()
            else:
                self.current_base += 1

    def _at_end(self):
        return self.current_base >= len(self.sequence)

    def _test_for_site_start(self):
        return self.p_not_bg[self.current_base] > self.threshold

    def _calculate_site(self):
        print 'Starting site at base %d' % self.current_base
        site_states = list()
        bases = list()
        last_state = self.ab[self.current_base].argmax()
        # find next state
        self.current_base += 1
        while not self._at_end() and last_state not in self.traits.background_states:
            # get prob of being in each state
            site_states.append(last_state)
            bases.append(self.sequence[self.current_base-1])
            p_next_state = self.beta[self.current_base] * self.A[last_state]
            last_state = p_next_state.argmax()
            self.current_base += 1
        return site_states, bases

# site_finder = SiteFinder(hmm, traits, threshold=0.9)
# sites = list(site_finder(numpy_seqs[1]))
# raise RuntimeError

def analyse_sequence_for_sites(sequence, hmm, traits, threshold=.9):
    LL, alpha, beta, c = hmm.forward_backward(sequence)
    ab = alpha * beta
    p_state = (ab.T / ab.sum(axis=1)).T
    p_not_bg = 1. - p_state[:,traits.background_states].sum(axis=1)
    best_base = p_not_bg.argmax()
    best_score = p_not_bg[best_base]
    start = best_base
    while start >= 0 and p_not_bg[start-1] > .5 * best_score:
        start -= 1
    end = best_base
    while end < len(sequence) and p_not_bg[end] > .5 * best_score:
        end += 1
    return best_score, sequence[start:end]


#
# Initialise the logging
#
logging.basicConfig(level=logging.INFO)



#
# Parse the options
#
option_parser = OptionParser()
option_parser.add_option(
  "-p",
  "--p-binding-site",
  dest="p_binding_site",
  default=0.01,
  type='float',
  help="Probability of a binding site in the model."
)
option_parser.add_option(
  "-m",
  "--models-file",
  dest="models_file",
  help="File in which the gapped PSSMs are stored."
)
option_parser.add_option(
  "-s",
  "--sequences-file",
  dest="sequences_file",
  help="FASTA file in which the sequences are stored."
)
option_parser.add_option(
  "--threshold-graph",
  dest="threshold_graph",
  help="file to write an image showing how # seqs with site varies by threshold."
)
# sys.argv='dummy.py -m /home/reid/T00759.pssm -s /home/reid/T00759.fa --threshold-graph test.png'.split()
options, args = option_parser.parse_args()
for option in option_parser.option_list:
    if option.dest:
        logging.info('%s: %s (%s)', option.dest, str(getattr(options, option.dest)), option.help)



logging.info('Loading sequences: %s', options.sequences_file)
sequences = dict(sequences_from_fasta(options.sequences_file))
numpy_seqs = dict((desc, seq_to_numpy(seq)) for desc, seq in sequences.iteritems())
logging.info('Loaded %d sequences', len(sequences))


logging.info('Parsing PSSMs: %s', options.models_file)
pssms = list(parse_models(open(options.models_file)))


logging.info('Building models')
models = [
  build_hmm_from_semi_parsed(parsed, p_binding_site=options.p_binding_site)
  for parsed in pssms
]

logging.info('Analysing sequences')
p_binding_sites = list()
sites_file = open('sites.txt', 'w')
for hmm, traits in models:
    for desc, sequence in numpy_seqs.iteritems():
        p_binding_site, site_seq = analyse_sequence_for_best_site(sequence, hmm, traits)
        p_binding_sites.append(p_binding_site)
        site = numpy_to_seq(site_seq)
        logging.info('%s: p(binding site)=%12g, sequence=%s', desc, p_binding_site, site)
        print >> sites_file, '%s, %12g, %s' % (desc, p_binding_site, site)
sites_file.close()

if options.threshold_graph:
    logging.info('Building threshold image: %s', options.threshold_graph)
    import pylab as P
    p_binding_sites.append(0.)
    p_binding_sites.sort()
    num_seqs = range(len(p_binding_sites))[::-1]
    P.figure()
    P.ylabel('# seqs with site')
    P.xlabel('p(binding) threshold')
    P.step(p_binding_sites, num_seqs, where='post')
    P.savefig(options.threshold_graph)
    P.close()
