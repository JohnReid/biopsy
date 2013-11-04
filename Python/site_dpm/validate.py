#
# Copyright John Reid 2009
#

"""
Code to validate transcriptional programs against certain groups of TFs and genes from
the literature.
"""

from shared import *
from generate_validation_sets import generate_validation_sets, get_kegg_sets, get_symatlas_sets, get_literature_sets
from tp_threshold import threshold_tps
from gene_set_enrichment import test_enrichment


def get_p_value_threshold_for(reference_name):
    return 1e-3

def q_value_analysis(p_values):
    logging.info('q-value analysis on %d p-values', len(p_values))
    from rpy2.robjects import r, FloatVector
    r.library('qvalue')
    ps = FloatVector(p_values)
    qobj = r.qvalue(ps)
    return qobj

def validate_tp_set(k, t, tp_set, reference_name, reference_set, universe, latex_f=None, p_values=None):
    test_tp_drawn, num_in_ref, num_not_in_ref, num_in_tp, p_value = test_enrichment(tp_set, reference_set, universe)
    assert len(tp_set) == num_in_tp
    assert len(reference_set) == num_in_ref
    if None != p_values:
        p_values.append(p_value)
    if p_value < get_p_value_threshold_for(reference_name):
        logging.info(
            'TP %3d; %s; % 4d/%4d in % 5d/%5d; %e; %s',
            k, t, test_tp_drawn, num_in_tp, num_in_ref, num_in_ref+num_not_in_ref, p_value, reference_name
        )
        if latex_f:
            print >> latex_f, '%3d & %s & %5d & %80s & %5d & %-5d & %.1e' % (
                k,
                t,
                num_in_tp,
                reference_name,
                test_tp_drawn,
                num_in_ref,
                p_value
            )


def validate_tp(tp, factor_validation_sets, target_validation_sets, factor_universe, target_universe):
    logging.info('Validating transcriptional program %d against %d sets of targets.', tp.k, len(target_validation_sets))


def restrict_validation_sets(validation_sets, universe):
    for vs in validation_sets:
        vs.intersection_update(universe)


@log_exceptions()
@caching_decorator('validate')
def validate():
    """
    Validate the transcriptional programs against the validation sets.
    """
    factor_validation_sets, target_validation_sets = generate_validation_sets()
    transcriptional_programs, factor_universe, target_universe = threshold_tps()
    factor_universe = set(factor_universe)
    target_universe = set(target_universe)
    restrict_validation_sets(factor_validation_sets.values(), factor_universe)
    restrict_validation_sets(target_validation_sets.values(), target_universe)
    latex_f = open(os.path.join(options.output_dir, 'validation.tex'), 'w')
    for tp in transcriptional_programs:
        #logging.info('Validating transcriptional program %d against %d sets of factors.', tp.k, len(factor_validation_sets))
        for name, reference_set in factor_validation_sets.iteritems():
            validate_tp_set(tp.k, 'Factors', tp.factors, name, reference_set, factor_universe, latex_f)
        #logging.info('Validating transcriptional program %d against %d sets of targets.', tp.k, len(target_validation_sets))
        for name, reference_set in target_validation_sets.iteritems():
            validate_tp_set(tp.k, 'Targets', tp.targets, name, reference_set, target_universe, latex_f)
    latex_f.close()


def validate_one_set(validation_sets, latex_f=None):
    transcriptional_programs, factor_universe, target_universe = threshold_tps()
    target_universe = set(target_universe)
    restrict_validation_sets(validation_sets.values(), target_universe)
    overall_p_values = list()
    for tp in transcriptional_programs:
        #logging.info('Validating transcriptional program %d against %d sets of targets.', tp.k, len(target_validation_sets))
        tp_p_values = list()
        for name, reference_set in validation_sets.iteritems():
            validate_tp_set(tp.k, 'Targets', tp.targets, name, reference_set, target_universe, latex_f, tp_p_values)
        overall_p_values.append(lou_jost_multiple_p_value_adjustment(reduce(float.__mul__, tp_p_values), len(tp_p_values)))
    return overall_p_values


def validate_one_factor_set(validation_sets, latex_f=None):
    transcriptional_programs, factor_universe, target_universe = threshold_tps()
    factor_universe = set(factor_universe)
    restrict_validation_sets(validation_sets.values(), factor_universe)
    overall_p_values = list()
    for tp in transcriptional_programs:
        #logging.info('Validating transcriptional program %d against %d sets of targets.', tp.k, len(target_validation_sets))
        tp_p_values = list()
        for name, reference_set in validation_sets.iteritems():
            validate_tp_set(tp.k, 'Factors', tp.factors, name, reference_set, factor_universe, latex_f, tp_p_values)
        overall_p_values.append(lou_jost_multiple_p_value_adjustment(reduce(float.__mul__, tp_p_values), len(tp_p_values)))
    return overall_p_values

def validiate_program_28():
    factor_validation_sets, target_validation_sets = generate_validation_sets()
    transcriptional_programs, factor_universe, target_universe = threshold_tps()
    factor_universe = set(factor_universe)
    target_universe = set(target_universe)
    restrict_validation_sets(factor_validation_sets.values(), factor_universe)
    restrict_validation_sets(target_validation_sets.values(), target_universe)
    tp = transcriptional_programs[28]
    for name, reference_set in target_validation_sets.iteritems():
        validate_tp_set(tp.k, 'Factors', tp.targets, name, reference_set, target_universe)

def lou_jost_multiple_p_value_adjustment(k, n):
    """
    U{http://www.loujost.com/Statistics%20and%20Physics/Significance%20Levels/CombiningPValues.htm}

    @arg k: product of all p-values
    @arg n: number of p-values
    """
    import math
    if n < 2:
        raise ValueError('Only works for n >= 2.')
    s = 0.
    minus_log_k = -math.log(k)
    for i in xrange(n-1):
        if 0 == i:
            minus_log_k_term = 1.
            fact = 1
        else:
            fact *= i
            minus_log_k_term *= minus_log_k
        s += minus_log_k_term / fact
    return s * k



if '__main__' == __name__:
    #validate()

    validiate_program_28()

    if False:
        latex_f = open(os.path.join(options.output_dir, 'validation-kegg.tex'), 'w')
        validation_sets = get_literature_sets()
        keys = set(validation_sets.keys())
        for k in keys:
            if not k.startswith('tremor'):
                del validation_sets[k]
        p_values = validate_one_set(validation_sets)
        q_obj = q_value_analysis(p_values)
        latex_f.close()
