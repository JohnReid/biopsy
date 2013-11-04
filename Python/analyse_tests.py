#
# Copyright John Reid 2006
#

import biopsy, time, pickle

num_thresholds = 1000

def only_matrices( pssm_list ):
    result = biopsy.SequenceVec()
    result.extend( [ p for p in pssm_list if p.startswith('M') ] )
    return result

pssms = biopsy.get_transfac_pssm_accessions(
        biopsy.get_default_transfac_pssm_filter() )
pssms = only_matrices( pssms )
print '# pssms:', len( pssms )

site_tcs = biopsy.get_transfac_site_test_cases()
fragment_tcs = biopsy.get_transfac_fragment_test_cases()
for tcs in [ site_tcs, fragment_tcs ]:
    for tc in tcs:
        tc.positives = only_matrices( tc.positives )
        tc.add_negatives_from( pssms )

def compare_params(
        classifiers,
        parameter_settings,
        test_cases,
        roc_points_filename
):
    try:
        roc_points = pickle.load( open(roc_points_filename, 'r') )
    except:
        roc_points = { }
    roc_points = { }

    for param_key in parameter_settings:
        print 'Params:', param_key
        for classifier_key, classifier in classifiers:
            roc_points[
                    classifier_key + '; ' + param_key
            ] = biopsy.evaluate_test_case_classifier(
                    test_cases,
                    classifier,
                    biopsy.test_case_truth,
                    num_thresholds = num_thresholds )

        pickle.dump(
                roc_points,
                open(roc_points_filename, 'w')
        )

    return roc_points


def show_roc_points( roc_points, title = None ):
    import pylab
    pylab.figure()
    for k, rp in roc_points.iteritems():
        biopsy.plot_roc_points( rp, label = str(k) )
    pylab.plot( [0,1], [0,1], 'k--', label = 'random' ) # straight line for random classifier
    pylab.legend( loc = 4 )
    if title: pylab.title( title )
    print 'Displaying ROC curves'
    pylab.show()

def compare_pseudo_counts():
    def gen_pseudo_counts():
        for pseudo_count in [
                0.00,
                0.01,
                0.10,
                0.25,
                0.50,
                0.75,
                1.00,
        ]:
            biopsy.PssmParameters.singleton().pseudo_counts = pseudo_count
            biopsy.clear_pssm_cache()
            yield "Pseudo-count: %.2f" % pseudo_count

    classifiers = [
            ( 'Bayesian', biopsy.NonPhyloClassifier( 0.0 ) ),
            ( 'Transfac min fp', biopsy.TransfacClassifier( biopsy.get_transfac_pssm_min_fp_threshold ) ),
    ]

    return compare_params(
            classifiers,
            gen_pseudo_counts(),
            fragment_tcs,
            'pseudo_counts_roc_points.pickle' )

def compare_cumulative():
    def gen_cumulative():
        for cumulative in [ True, False ]:
            for p_value in [ True, False ]:
                biopsy.PssmParameters.singleton().use_cumulative_dists = cumulative
                biopsy.PssmParameters.singleton().use_p_value = p_value
                yield "%s - %s" % (
                        cumulative and 'Cumulative' or 'Exact',
                        p_value and 'p-value' or 'bayesian',
                )

    classifiers = [
            ( ' - ', biopsy.NonPhyloClassifier( 0.0 ) ),
    ]

    return compare_params(
            classifiers,
            gen_cumulative(),
            fragment_tcs,
            'cumulative_roc_points.pickle' )

pseudo_count_roc_points = compare_pseudo_counts()
cumulative_roc_points = compare_cumulative()

show_roc_points( pseudo_count_roc_points )
show_roc_points( cumulative_roc_points )
