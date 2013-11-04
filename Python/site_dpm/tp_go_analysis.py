#
# Copyright John Reid 2009
#

"""
Code to run GO enrichment analysis on transcriptional programs.
"""

from shared import *
import tp_threshold, go, sys, topgo

@log_exceptions()
@caching_decorator('go-analysis')
def go_analysis():
    logging.info('Running GO analysis: p-value threshold=%e; topGO method=%s', options.go_p_value_threshold, options.topgo_method)
    transcriptional_programs, factor_universe, target_universe = tp_threshold.threshold_tps()
    genes_2_GO, go_context = go.initialise_go_context(factor_universe, target_universe, options.go_ontologies)
    go_analyses = list()
    f = open(os.path.join(options.output_dir, 'go-analyses.txt'), 'w')
    for tp in transcriptional_programs:
        go_analysis = go.TPGoAnalysis(tp, go_context, options.go_p_value_threshold, options.topgo_method)
        go_analysis.print_go_analyses(f, options.go_p_value_threshold, log=True)
        go_analysis.print_go_analyses(sys.stdout, options.go_p_value_threshold)
        go_analyses.append(go_analysis)
    f.close()
    write_latex(go_analyses)


def write_latex(go_analyses):
    latex_f = open(os.path.join(options.output_dir, 'go-analyses.tex'), 'w')
    print_latex(go_analyses, f=latex_f, factors=True)
    print_latex(go_analyses, f=latex_f, factors=False)
    latex_f.close()



def print_latex(go_analyses, f=sys.stdout, factors=True):
    transcriptional_programs, factor_universe, target_universe = tp_threshold.threshold_tps()
    t = factors and 'Factors' or 'Targets'
    print >> f, '%% %s' % t
    print >> f, 'TP & %s & GO term & & GO description & \\multicolumn{2}{c}{annotated} & $p$-score \\\\' % t
    print >> f, '\\hline'
    for k, (tp, tp_analysis) in enumerate(zip(transcriptional_programs, go_analyses)):
        if factors:
            type_analysis = tp_analysis.factors_go_analysis
            size = len(tp.factors)
        else:
            type_analysis = tp_analysis.targets_go_analysis
            size = len(tp.targets)
        for ontology, analysis in type_analysis.iteritems():
            if None != analysis:
                for go_id, go_term, annotated, significant, expected, pvalue in topgo.yield_stats(analysis):
                    # Program & Factors & GO term & & GO description & \multicolumn{2}{c}{annotated} & $p$-score \\
                    print >> f, '% 4d & %5d & %s & %s & %42s & % 4d & %-4d & %.1e \\\\' % (
                        tp.k,
                        size,
                        go_id,
                        ontology,
                        go_term,
                        significant,
                        annotated,
                        pvalue
                    )
    print >> f, '\\\\'

if '__main__' == __name__:
    go_analysis()
