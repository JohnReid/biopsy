#
# Copyright John Reid 2009
#


"""
Do a bootstrap analysis of the GO p-values.
"""


from shared import *
import tp_threshold, go, topgo, numpy, os
from infpy.bootstrap import *
from cookbook.pylab_utils import pylab_ioff


@log_exceptions()
@caching_decorator('bootstrap')
def do_bootstrap():
    """
    Do a bootstrap analysis on GO p-values.
    """
    logging.info('Running GO bootstrap analysis with %d samples: topGO method=%s', options.num_bootstrap_samples, options.topgo_method)
    transcriptional_programs, factor_universe, target_universe = tp_threshold.threshold_tps()
    genes_2_GO, go_context = go.initialise_go_context(factor_universe, target_universe, options.go_ontologies)
    tp_sizes = filter(None, map(len, (tp.targets for tp in transcriptional_programs)))
    p_values = list()
    for sample in generate_bootstrap_samples(options.num_bootstrap_samples, target_universe, tp_sizes):
        go_analysis = dict(
            (
                ontology,
                go.try_go_analysis(
                    go_data,
                    sample,
                    1.,
                    options.topgo_method
                )
            )
            for ontology, go_data
            in go_context.targets_go_data.iteritems()
        )
        best_p_value = min(map(topgo.p_value_from_r, map(topgo.best_p_value, go_analysis.values())))
        p_values.append(best_p_value)
    logging.info('GO bootstrap analysis completed')
    return p_values




@pylab_ioff
def plot_bootstrap():
    import pylab as P
    p_values = do_bootstrap()
    bootstrap_dir = os.path.join(options.output_dir, 'bootstrap')
    ensure_dir_exists(bootstrap_dir)

    # p-value histogram
    fig = P.figure()
    P.hist(p_values, bins=24)
    P.title('Bootstrap p-values')
    P.xlabel('p-values')
    P.savefig(os.path.join(bootstrap_dir, 'p-values.png'))
    P.savefig(os.path.join(bootstrap_dir, 'p-values.eps'))
    P.close(fig)

    # log p-values
    fig = P.figure()
    P.title('Bootstrap p-values')
    P.xlabel('log p-values')
    P.hist(numpy.log10(p_values), bins=24)
    P.savefig(os.path.join(bootstrap_dir, 'log-p-values.png'))
    P.savefig(os.path.join(bootstrap_dir, 'log-p-values.eps'))
    P.close(fig)



if '__main__' == __name__:
    #do_bootstrap()
    plot_bootstrap()
