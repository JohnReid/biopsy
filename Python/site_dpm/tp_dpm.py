#
# Copyright John Reid 2009
#

"""
Code to deal with transcriptional program DPMs.
"""

from shared import *
from hdpm import HDPM, DpmSummariser, DpmInferenceHistory
from infpy.convergence_test import LlConvergenceTest
from cookbook.pylab_utils import pylab_ioff


def create_dpm(analysis, K):
    """
    Takes the output from the remome analysis and builds the Dirichlet process mixture.
    """
    from itertools import repeat
    #
    # Work out which factors we have across all our remos so we can map factors to indices
    #
    factors = set()
    for hit_counts in analysis.values():
        factors.update(hit_counts.keys())
    factors = list(factors) # make into a list so we can index it
    # map each factor to an index
    factor_indices = dict((factor, i) for i, factor in enumerate(factors))
    logging.info('Found %d distinct factors in binding sites', len(factors))

    targets = list(k for k in analysis.keys() if len(analysis[k]))
    target_indices = dict((target, i) for i, target in enumerate(targets))

    def hit_counts_to_document(hit_counts):
        result = list()
        for factor, count in hit_counts.iteritems():
            result.extend(repeat(factor_indices[factor], count))
        return result

    dpm_input = [
      numpy.array(hit_counts_to_document(analysis[target]))
      for target
      in targets
    ]

    logging.info('a_alpha: %3g; b_alpha: %3g', options.a_alpha, options.b_alpha)
    logging.info('a_beta:  %3g; b_beta:  %3g', options.a_beta, options.b_beta)
    logging.info('a_gamma: %3g; b_gamma: %3g', options.a_gamma, options.b_gamma)

    logging.info('Creating DPM with %d transcriptional programs', K)
    dpm = HDPM(
      documents=dpm_input,
      K=K,
      W=len(factors),
      a_alpha=options.a_alpha,
      b_alpha=options.b_alpha,
      a_beta=options.a_beta,
      b_beta=options.b_beta,
      a_gamma=options.a_gamma,
      b_gamma=options.b_gamma,
      a_tau=options.a_tau
    )

    ensembl_names = get_all_ensembl_names()
    target_names = [ensembl_names.get(str(g), '<unknown>') for g in targets]
    factor_names = [ensembl_names[str(g)] for g in factors]
    summariser = DpmSummariser(
      dpm,
      os.path.join(summaries_dir, 'dpm-summary'),
      target_names,
      factor_names,
      document_tag='target',
      topic_tag='program',
      word_tag='factor',
      occurence_tag='binding event'
    )

    return factors, factor_indices, targets, target_indices, dpm_input, dpm, summariser, ensembl_names



def infer_dpm(dpm, summariser, min_iters=10, max_iters=2):
    #
    # update DPM to convergence
    #
    convergence_test = LlConvergenceTest(should_increase=False, use_absolute_difference=True)
    start = time()
    summariser.log_static_info()
    summariser.log_dynamic_info()
    logging.info('Updating Dirichlet process mixture until log likelihood converges or %d iterations', max_iters)
    history = DpmInferenceHistory(dpm, summariser)
    for i in xrange(max_iters):
        dpm.update()
        summariser.log_dynamic_info()
        LL = history.iteration()
        logging.info('iteration %3d; LL: %f', i, LL)
        if convergence_test(LL) and i >= min_iters:
            break
    history.plot()
    total_elapsed = time()-start
    logging.info('%d iterations took %f secs, %f secs/iteration', i + 1, total_elapsed, total_elapsed/(i+1))
    logging.info('LL improved from %f to %f', convergence_test.LLs[0], convergence_test.LLs[-1])
    return convergence_test, history






def output_summary(summariser):
    summariser.plot_num_documents_by_top_topics()
    summariser.topic_sizes()
    summariser.topic_sizes_log_scale()
    summariser.histograms()
    summariser.make_heat_maps()
    summariser.plot_word_document_scatter()



def write_background_sets(output_dir, factors, targets):
    """
    Write background gene/pssm sets
    """
    bg_pssm_file = os.path.join(output_dir, 'background_factors.txt')
    logging.info('Writing background factors to: %s', bg_pssm_file)
    open(bg_pssm_file, 'w').write('\n'.join(str(g) for g in factors))

    bg_targets_file = os.path.join(output_dir, 'background_targets.txt')
    logging.info('Writing background targets to: %s', bg_targets_file)
    open(bg_targets_file, 'w').write('\n'.join(str(g) for g in targets))


@pylab_ioff
def plot_hits_per_target(hit_counts):
    import pylab as P
    P.figure()
    plot_hit_counts_histogram(hit_counts, bins=100)
    P.savefig(os.path.join(summaries_dir, 'num-hits-per-target.png'))
    P.savefig(os.path.join(summaries_dir, 'num-hits-per-target.eps'))
    P.close()


@pylab_ioff
def plot_hit_counts_histogram(hit_counts, *args, **kwds):
    import pylab as P
    P.hist(num_hits_per_target(hit_counts).values(), *args, **kwds)
    P.title('# hits per target')
    P.xlabel('# hits')
    P.ylabel('# targets')


@pylab_ioff
def plot_programs_info(summariser):
    for k in xrange(summariser.statistics.num_topics_used):
        import pylab as P
        fig = summariser.plot_factor_enrichment(k)
        P.savefig(os.path.join(get_programs_dir(), '%03d-factor-enrichment.png' % k))
        P.close(fig)


def num_hits_per_target(hit_counts):
    "@return: A dict mapping targets to numbers of hits."
    return dict((g, sum(hits.values())) for g, hits in hit_counts.iteritems())


@log_exceptions()
@caching_decorator('dpm')
def create_and_infer_dpm():
    #hit_counts = get_hit_counts(
        #remome_threshold=remome_threshold,
        #masked=masked,
        #use_max_chain=use_max_chain,
        #analysis_threshold=analysis_threshold,
    #)
    import hit_counts
    hit_counts = hit_counts.get_hit_counts()
    num_binding_sites = sum(sum(seq_counts.values()) for seq_counts in hit_counts.values())
    logging.info('Have %d putative binding sites spread across %d targets', num_binding_sites, len(hit_counts))
    plot_hits_per_target(hit_counts)
    factors, factor_indices, targets, target_indices, dpm_input, dpm, summariser, ensembl_names = create_dpm(hit_counts, options.K)
    write_background_sets(options.output_dir, factors, targets)
    convergence_test, history = infer_dpm(dpm, summariser, min_iters=options.min_iters, max_iters=options.max_iters)
    summariser.statistics.document_topic_threshold = options.document_topic_threshold
    summariser.statistics.topic_word_threshold = options.topic_word_threshold
    summariser.statistics.update()
    output_summary(summariser)
    plot_programs_info(summariser)
    return factors, factor_indices, targets, target_indices, dpm_input, dpm, summariser, ensembl_names

if '__main__' == __name__:
    factors, factor_indices, targets, target_indices, dpm_input, dpm, summariser, ensembl_names = create_and_infer_dpm()
