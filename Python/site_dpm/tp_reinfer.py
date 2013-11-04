#
# Copyright John Reid 2009
#

"""
Takes an already inferred DPM and does more iterations.
"""


from shared import *


@log_exceptions()
@caching_decorator('re-inferred-dpm')
def reinfer_dpm():
    import tp_dpm
    factors, factor_indices, targets, target_indices, dpm_input, dpm, summariser, ensembl_names = tp_dpm.create_and_infer_dpm()
    convergence_test, history = tp_dpm.infer_dpm(dpm, summariser, min_iters=options.min_iters, max_iters=options.max_iters)
    tp_dpm.output_summary(summariser)
    tp_dpm.plot_programs_info(summariser)
    return factors, factor_indices, targets, target_indices, dpm_input, dpm, summariser, ensembl_names


if '__main__' == __name__:
    reinfer_dpm()
