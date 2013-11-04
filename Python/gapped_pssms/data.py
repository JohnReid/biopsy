#
# Copyright John Reid 2008
#

"""
Code to locate and retrieve data for gapped PSSM work.
"""

import os, os.path

all_fragments = [
  'T00594',
  'T00163',
  'T00759',
  'T00368',
  'T03286',
  'T00140',
  'T00671',
  'T09363',
  'T03828',
  'T08969',
  'T00781',
]
"""
The chip-chip fragment names arranged in order of size.
"""

test_set_fragments = [
  'T00594',
  'T00163',
  'T00368',
  'T03286',
  'T00140',
  'T00671',
  'T00759',
]
"""
The chip-chip fragment names that have test sets associated with them arranged in order of size.
"""

cross_folds = 5

def _looks_like_data_dir(candidate_dir):
    "@return: True if and only if candidate_dir looks like the data directory."
    return os.access(os.path.join(candidate_dir, 'T00594.fa'), os.R_OK)

def _try_candidate_dirs(candidate_dirs):
    "@return: The first candidate directory in the sequence that looks like the data directory."
    for candidate_dir in candidate_dirs:
        if _looks_like_data_dir(candidate_dir):
            return candidate_dir

if 'nt' == os.name:
    candidate_dirs = [
      os.path.join('C:\\', 'Analysis', 'GappedPssms', 'fragments')
    ]
else:
    candidate_dirs = [
      '/home/reid/Analysis/GappedPssms/fragments/fasta',
      '/home/john/Analysis/GappedPssms/fragments/fasta',
    ]

data_dir = _try_candidate_dirs(candidate_dirs)
if None == data_dir:
    import warnings
    warnings.warn('Could not locate data for gapped PSSM work.')

if None != data_dir:
    results_dir = os.path.join(data_dir, 'results')
    test_sets_dir = os.path.join(data_dir, 'test-sets')
    cross_validate_dir = os.path.join(data_dir, 'revised-fragments', 'fasta')
else:
    results_dir = None
    test_sets_dir = None
    cross_validate_dir = None



def fasta_file_for_fragment(fragment):
    "@return: The FASTA filename for the fragment data."
    return os.path.join(data_dir, '%s.fa' % fragment)

def fasta_files_for_fragment_cross_fold(fragment, cross_fold_idx):
    "@return: The FASTA filenames for the training and test fragment data."
    fragment_dir = os.path.join(test_sets_dir, fragment)
    return (
      os.path.join(cross_validate_dir, '%strim-train-x%d.fa' % (fragment, cross_fold_idx)),
      os.path.join(cross_validate_dir, '%strim-test-x%d.fa' % (fragment, cross_fold_idx)),

    )

def fasta_file_for_synthetic_data(tag):
    "@return: The FASTA filename for the sythetic data."
    return os.path.join(data_dir, 'synthetic', 'synthetic-sequences-%s.fa' % tag)

def training_test_sequences():
    """
    @return: A list of (training, test) pairs of sequences
    """
    from gapped_pssms import sequence
    sequences = []
    for fragment in test_set_fragments:
        for cross_fold in xrange(cross_folds):
            sequences.append(
              map(
                sequence.convert_fasta_sequences,
                fasta_files_for_fragment_cross_fold(fragment, cross_fold)
              )
            )
    return sequences
