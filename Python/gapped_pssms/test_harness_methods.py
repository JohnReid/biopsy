#
# Copyright John Reid 2009
#


"""
Code to examine directories for the output of methods.
"""


import os, logging, sys, glam2, re, glob
from optparse import OptionParser
from test_harness_2 import TestHarness
from itertools import imap

#
# Initialise the logging
#
format='%(asctime)-15s %(message)s'
logging.basicConfig(level=logging.INFO, format=format)
logging.info('Command line: %s', ' '.join(sys.argv))
logging.info('Current working directory: %s', os.getcwd())



class Glam2Dir(object):

    def __init__(self, dir):
        self.dir = dir
        self.files = dict()
        for best in glam2.best_files(dir):
            fragment, cross_fold_index, seed = glam2.interpret_output_filename(best)
            self.files[(fragment, cross_fold_index)] = best


class REDir(object):

    @staticmethod
    def interpret_output_filename(re, filename):
        "Parse the information in the filename."
        match = re.search(filename)
        if not match:
            return None
        fragment, cross_fold_index = match.groups()
        cross_fold_index = int(cross_fold_index)
        return fragment, cross_fold_index


    def __init__(self, dir, re):
        self.dir = dir
        self.files = dict()
        for filename in glob.glob(os.path.join(dir, '*.pssm')):
            key = REDir.interpret_output_filename(re, filename)
            if key:
                self.files[key] = filename


class MemeDir(REDir):
    re = re.compile('vm-(.*)-motif-h[0-9]+-v[0-9]+-x([0-9]+).pssm')
    def __init__(self, dir):
        REDir.__init__(self, dir, MemeDir.re)


class GappedDir(REDir):
    re = re.compile('([^/]*)-([0-9]+).pssm')
    def __init__(self, dir):
        REDir.__init__(self, dir, GappedDir.re)


def load_dir(dir):
    "Examines directory for method output and loads it."
    types = [
        MemeDir,
        Glam2Dir,
        GappedDir,
    ]
    for t in types:
        d = t(dir)
        if len(d.files):
            return d
    raise RuntimeError('Could not find method output in %s' % dir)



if '__main__' == __name__:
    analysis_dir = '/home/john/Analysis/GappedPssms'
    glam2_dir = load_dir(os.path.join(analysis_dir, 'GLAM2/cross-validate'))
    print glam2_dir.files
    meme_dir = load_dir(os.path.join(analysis_dir, 'MEME/x-validate'))
    print meme_dir.files
    gapped_dir = load_dir(os.path.join(analysis_dir, 'Gapped'))
    print gapped_dir.files
