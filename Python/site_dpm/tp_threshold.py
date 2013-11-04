#
# Copyright John Reid 2009
#

"""
Code to threshold transcriptional programs.
"""

from shared import *
import tp_dpm
import basic_tp


@log_exceptions()
@caching_decorator('tps')
def threshold_tps():
    factor_universe, factor_indices, target_universe, target_indices, dpm_input, dpm, summariser, ensembl_names = tp_dpm.create_and_infer_dpm()
    transcriptional_programs = [
        basic_tp.tp_from_dpm_summary(summariser, factor_universe, target_universe, k)
        for k in xrange(summariser.statistics.num_topics_used)
    ]
    for tp in transcriptional_programs:
        tp.write_files(ensembl_names)
    return transcriptional_programs, factor_universe, target_universe


if '__main__' == __name__:
    threshold_tps()
