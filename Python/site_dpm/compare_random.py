#
# Copyright John Reid 2008, 2009
#

"""
Code to compare random GO analyses.
"""

from shared import *


def histogram_p_values(p_values, tag):
    import pylab as P
    interactive_mode = P.isinteractive()
    P.ioff()
    try:
        P.figure()
        P.hist(p_values)
        P.title('%s p-values for transcriptional programs' % tag)
        P.xlabel('p-values')
        P.savefig(os.path.join(output_dir, 'p-values-for-%ss.png' % tag), format='PNG')
        P.close()
    finally:
        if interactive_mode:
            P.ion()


def random_samples_like(transcriptional_programs):
    'Sample sets of targets and factors randomly that are the same size as those in the programs.'
    from rpy import r
    target_samples = []
    factor_samples = []
    for tp in transcriptional_programs:
        if tp.targets_go_analysis:
            target_samples.append(r.sample(tp.genes, len(tp.tp_targets)))
        if tp.factors_go_analysis:
            factor_samples.append(r.sample(tp.factors, len(tp.tp_factors)))
    return target_samples, factor_samples


def random_go_analyses(samples, go_data, p_value_threshold):
    return [GoAnalysis(go_data, sample, p_value_threshold) for sample in samples]

def print_analyses(samples, go_analyses, tag):
    for k, (sample, go_analysis) in enumerate(zip(samples, go_analyses)):
        print '************ Random sample: %d **************; # %s=%d' % (k, tag, len(sample))
        print str(go_analysis.results)

def best_p_value(go_analysis):
    return float(go_analysis.results_unthresholded.classic[0])

def calculate_best_p_value(go_data, sample):
    from rpy import RPyRException
    try:
        return float(GoAnalysis(go_data, sample, 0.01).results_unthresholded.classic[0])
    except RPyRException:
        return 1.0

def generate_best_p_values(sample_sizes, num_samples_at_each_size, universe, go_data):
    from rpy import r
    for size in sample_sizes:
        yield size, [
          log10(calculate_best_p_value(go_data, r.sample(universe, size)))
          for i in xrange(num_samples_at_each_size)
        ]

def create_p_value_boxplot_eps(best_p_values, filename):
    from rpy import r
    r.postscript(filename, horizontal=False, height=4.5, width=6, pointsize=10)
    try:
        keys = best_p_values.keys()
        keys.sort()
        r.boxplot(
            map(best_p_values.get, keys), 
            names=map(str, keys),
            xlab="sample size", ylab="p-score")
    finally:
        r.dev_off()

def generate_best_p_value_eps(sample_sizes, num_samples_at_each_size, universe, go_data, filename):
    best_p_values = dict(generate_best_p_values(sample_sizes, num_samples_at_each_size, universe, go_data))
    create_p_value_boxplot_eps(best_p_values, filename)
    return best_p_values


if '__main__' == __name__:
    1/0
    num_samples_at_each_size = 100
    #num_samples_at_each_size = 2
    factor_best_p_values = generate_best_p_value_eps(
      xrange(1,51),
      num_samples_at_each_size,
      factors,
      factors_go_data,
      os.path.join(output_dir, "factor_random_samples_best_p_values.eps")
    )
    target_best_p_values = generate_best_p_value_eps(
      xrange(2,202,4),
      num_samples_at_each_size,
      map(str, genes),
      targets_go_data,
      os.path.join(output_dir, "target_random_samples_best_p_values.eps")
    )
    cPickle.dump((factor_best_p_values, target_best_p_values), open('best-p-values.pickle', 'w'))
    #factor_best_p_values, target_best_p_values = cPickle.load(open('best-p-values.pickle'))

    # get a list of all the p-values
    all_random_p_values = list(chain(chain(*factor_best_p_values.values()), chain(*target_best_p_values.values())))
    all_random_p_values.sort()
    assert len(all_random_p_values) == 100 * (len(factor_best_p_values) + len(target_best_p_values))



    # fit extreme values
    def dlgumbel(x, mu, sigma):
        z = -(x-mu)/sigma
        return z - exp(z) - log(sigma)

    def log_likelihood_random_values(params):
        return dlgumbel(-x, *params).sum()

    def to_optimise(params):
        return -log_likelihood_random_values(params)

    x = numpy.asarray(all_random_p_values)
    start_params = numpy.array([ 1., 1. ])
    import scipy.optimize as O
    xopt, fopt, func_calls, grad_calls, warnflag = O.fmin_cg(to_optimise, start_params, full_output=True)

    # we don't need what's below any more
    raise RuntimeError()




    L = len(all_random_p_values)
    proportions=[math.log10((i+1)/float(L)) for i, p_value in enumerate(all_random_p_values)]
    figure()
    plot(proportions, all_random_p_values)
    ylabel("p-score")
    xlabel("log(proportion of better p-scores)")
    savefig(os.path.join(output_dir, "random-p-score-proportions.eps"), format="EPS")
    savefig(os.path.join(output_dir, "random-p-score-proportions.png"), format="PNG")

    # plot not on log-scale - not as convincing
    #figure()
    #plot(numpy.exp(proportions), numpy.exp(all_random_p_values))




    factor_sizes = [len(tp.tp_factors) for tp in transcriptional_programs]
    target_sizes = [len(tp.tp_targets) for tp in transcriptional_programs]


    # Get the go analyses from the transcriptional programs
    factor_go_analyses = [tp.factors_go_analysis for tp in transcriptional_programs if tp.factors_go_analysis]
    factor_best_p_values = map(math.log10, imap(best_p_value, factor_go_analyses))

    # get best p-value from each analysis of programs
    target_go_analyses = [tp.targets_go_analysis for tp in transcriptional_programs if tp.targets_go_analysis]
    target_best_p_values = map(math.log10, imap(best_p_value, target_go_analyses))

    # Get random samples of same size as programs
    target_samples, factor_samples = random_samples_like(transcriptional_programs)

    # analyse random samples for GO term enrichment
    random_target_go_analyses = random_go_analyses(target_samples, targets_go_data, p_value_threshold)
    random_factor_go_analyses = random_go_analyses(factor_samples, factors_go_data, p_value_threshold)

    # take best p-value from each random analysis
    random_target_best_p_values = map(math.log10, imap(best_p_value, random_target_go_analyses))
    random_factor_best_p_values = map(math.log10, imap(best_p_value, random_factor_go_analyses))

    # scatter plot the 2 sets of p-values against each other and a straight line
    P.figure()
    P.scatter(target_best_p_values, random_target_best_p_values, c='blue')
    P.scatter(factor_best_p_values, random_factor_best_p_values, c='yellow')
    range = (max(gca().get_xlim()[0], gca().get_ylim()[0]), min(gca().get_xlim()[1], gca().get_ylim()[1]))
    P.plot(range, range, 'k:')
    P.xlabel('p-values from programs')
    P.ylabel('p-values from random samples')
    summariser.save_fig('p-value-compare-random-scatter')

    # cPickle.dump((target_best_p_values, factor_best_p_values, random_target_best_p_values, random_factor_best_p_values), open('best-p-value-analysis.pickle', 'w'))
