#
# Copyright John Reid 2009, 2010, 2011
#

"""
Default options for transcriptional program work.
"""


r_console_width = 240

# Remome options
remome_threshold = 100
masked = True

# UCSC options
use_ucsc_seqs = True
ucsc_use_masked_seqs = True

# Analysis options
use_max_chain=False
analysis_threshold = 0.05

# hit count options
hit_count_threshold = .25

# DPM options
min_iters=10
max_iters=120
a_alpha = 10000.
b_alpha = 1000.
a_beta = 10000.
b_beta = 10000.
a_gamma = 10000.
b_gamma = 10000.
a_tau = None
K = 600

# Posterior options
document_topic_threshold = topic_word_threshold = 5.

# output directory
output_dir = '.'

# SymAtlas options
symatlas_fold_change_threshold = 10.

# GO analysis options
go_p_value_threshold = 1e-5
go_ontologies = ('BP', 'MF', 'CC')
#go_ontologies = ('BP',)
go_evidence_codes_to_ignore=set()
topgo_method = 'classic'
num_bootstrap_samples = 100
