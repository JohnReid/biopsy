#
# Copyright John Reid 2008
#

"""
Code to rank PSSMs by "interesting-ness".

Information content.
Low-order predictability.
Number of sequences with sites.
"""

from gapped_pssms.parse_gapped_format import parse_models
from itertools import imap
from gapped_pssms.pssm_score import *
from cookbook import DictOfLists
import glob, logging, shutil, os, re
import hmm.pssm.logo as logo


def calculate_emissions(model):
    emissions = numpy.zeros((model.N, model.M))
    for i in xrange(model.N):
        assert model.emissions[i][0] == i
        emissions[i] = model.emissions[i][1]
    return emissions

def calculate_gap_probs(model):
    gap_probs = numpy.ones((model.N))
    for f, t, p in model.transitions:
        gap_probs[t] = p
    return gap_probs

class Log(object):
    """
    Parses log files.
    """

    log_file_name_re = re.compile('(.*).log')
    pssm_num_re = re.compile('PSSM ([0-9]+)')
    sites_re = re.compile('[Ff]ound ([0-9]+) sites. ([0-9]+)/([0-9]+) sequences have at least one site')

    def __init__(self, log_file):
        """
        Parses log files:

        ************** PSSM 4 **************
        Seed ctgctgtg with gap at 3 had 79 hits in 72/601 sequences
        Seed score: 2528.810289
        Found 238 sites. 145/601 sequences have at least one site
        Entropy/base        : 0.923442
        Information content : 10.238500
        """
        logging.info('Parsing log file %s', log_file)
        self.log_file = log_file
        self.site_numbers = dict()

        re_match = Log.log_file_name_re.match(os.path.basename(log_file))
        self.tag = re_match.group(1)
        logging.info('%s: %s', self.log_file, self.tag)

        for line in open(log_file):
            m = Log.pssm_num_re.search(line)
            if m:
                pssm_num = int(m.group(1))
                # logging.info('PSSM: %d', pssm_num)
            m = Log.sites_re.search(line)
            if m and -1 == line.find('Trained model'):
                num_sites = int(m.group(1))
                num_seqs_with_site = int(m.group(2))
                num_seqs = int(m.group(3))
                # logging.info('# sites: %d; # seqs with sites: %d; # seqs: %d', num_sites, num_seqs_with_site, num_seqs)
                self.site_numbers[pssm_num] = (num_sites, num_seqs_with_site, num_seqs)



class Pssm(object):
    pssm_file_name_re = re.compile('(.*)-([0-9]+).pssm')

    def __init__(self, pssm_file, log):
        self.pssm_file = pssm_file
        self.png_file = pssm_file.replace('.pssm', '.png')
        self.eps_file = pssm_file.replace('.pssm', '.eps')

        re_match = Pssm.pssm_file_name_re.match(os.path.basename(pssm_file))
        self.tag = re_match.group(1)
        self.pssm_idx = int(re_match.group(2))
        self.num_sites, self.num_seqs_with_site, self.num_seqs = log.site_numbers[self.pssm_idx]
        # logging.info('%s: %s %d %d', self.pssm_file, self.fragment, self.cross_fold, self.pssm_idx)

        self.model = parse_models(open(self.pssm_file)).next()
        self.emissions = calculate_emissions(self.model)
        self.gap_probs = calculate_gap_probs(self.model)
        self.first_order_entropy_score = calculate_first_order_entropy_score(self.emissions)
        self.information_content_score = calculate_information_content_score(self.emissions)
        self.num_seqs_with_site_score = float(self.num_seqs_with_site) / float(self.num_seqs)
        self.overall_score = weighted_geometric_mean(
          (self.first_order_entropy_score, self.information_content_score, self.num_seqs_with_site_score),
          [1.5                           , 1.                            , 1.]
        )
        logging.info(
          '%s; %8g; %8g; %8g; %8g',
          self.pssm_file,
          self.first_order_entropy_score,
          self.information_content_score,
          self.num_seqs_with_site_score,
          self.overall_score
        )

    def write_image(self):
        image = logo.pssm_as_image(
          self.emissions,
          transparencies=self.gap_probs
        )
        image.save(self.png_file, "PNG")
        image.save(self.eps_file, "EPS")



def montage(input_files, output_file):
    montage_cmd = 'montage -tile 1x -geometry x240 %s %s' % (' '.join(input_files), output_file)
    os.system(montage_cmd)


class PssmSet(object):
    def __init__(self, basename):
        self.basename = basename
        self.tag = os.path.basename(basename)
        self.log = Log('%s.log' % self.basename)
        self.pssms = dict(
          (num, Pssm('%s-%03d.pssm' % (self.basename, num), self.log))
          for num in self.log.site_numbers.keys()
        )

    def sorted_by_score(self):
        """
        Returns a list of pssms sorted by score.
        """
        sorted_pssms = self.pssms.values()
        sorted_pssms.sort(key=lambda p: p.overall_score, reverse=True)
        logging.info(' '.join(imap(str, (p.pssm_idx for p in sorted_pssms))))
        return sorted_pssms

    def montage_by_score(self):
        sorted_pssms = self.sorted_by_score()
        ranked_files = [p.png_file for p in sorted_pssms]
        ranked_file = '%s-ranked.png' % self.basename
        montage(ranked_files, ranked_file)


if '__main__' == __name__:
    logging.basicConfig(level=logging.DEBUG)

    import sys

    root_dir = sys.argv[1]
    tag_re = re.compile('(T.*).log')
    tags = map(lambda m: m.group(1), filter(None, imap(tag_re.search, glob.glob(os.path.join(root_dir, 'T*.log')))))
    # tags = ['T00140-3']
    for tag in tags:
        logging.info(tag)
        pssm_set = PssmSet(os.path.join(root_dir, tag))
        pssm_set.montage_by_score()
