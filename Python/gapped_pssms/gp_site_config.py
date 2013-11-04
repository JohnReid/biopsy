#
# Copyright John Reid 2009
#

import os

default_positive_sequence_pattern = '%strimRM.fa'
default_positive_sequence_train_pattern = '%strimRM-train-x%d.fa'
default_positive_sequence_test_pattern = '%strimRM-test-x%d.fa'

machine_name = os.uname()[1]
if 'john-dell' == machine_name:
    default_positive_sequence_dir = '/home/john/Analysis/GappedPssms/fragments/fasta/'
    default_background_sequence_dir = default_positive_sequence_dir

elif 'sysbio' == machine_name:
    default_positive_sequence_dir = '/home/reid/Analysis/GappedPssms/fragments/revised-fragments/fasta'
    default_background_sequence_dir = default_positive_sequence_dir

elif 'pc118.bio.warwick.ac.uk' == machine_name:
    default_positive_sequence_dir = '/home/kenneth/jan-2009'
    default_background_sequence_dir = '/home/kenneth/dec-2008'
