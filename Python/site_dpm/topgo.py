#
# Copyright John Reid 2009
#

"""
Code to interface to the topGO library.
"""

import rpy2.robjects, logging, os
from rpy2.robjects import r
from shared import *

r_str = r.str
r_print = r['print']


def set_r_width(width):
    "Set the R width option (useful for prettier printing)"
    r['options'](width=width)


def as_list(obj):
    if isinstance(obj, list):
        return obj
    return list(obj)

def go_enrichment_script_filename():
    "@return: The filename of the R GO enrichment script"
    try:
        return os.path.join(os.path.dirname(__file__), 'go-enrichment.R')
    except:
        return 'go-enrichment.R'


def initialise_topgo():
    "Initialise the library."
    r('source("%s")' % go_enrichment_script_filename())

def pretty_str(obj):
    "@return: A pretty printed string of the results."
    #a = with_mode(NO_CONVERSION, lambda: r.textConnection('tmpobj', 'w'))()
    conn = r.textConnection('tmpobj_pretty_str', 'w')
    r.sink(file=conn, type='output')
    r['print'](obj)
    r.sink()
    r.close(conn)
    return '\n'.join(r['tmpobj_pretty_str'])
    #str = with_mode(BASIC_CONVERSION, lambda: r('tmpobj'))()
    #return '\n'.join(as_list(str))

def convert_genes_2_GO(genes_2_GO):
    "@return: An R object that can be used in to map genes to GO identifiers."
    return r.list(**genes_2_GO)

def create_go_data(gene_names, r_genes_2_GO, ontology):
    "Create the topGO data to use to analyse transcriptional programs."
    logging.info('Creating topGO data for %d genes using ontology %s' % (len(gene_names), ontology))
    return r.buildGoData(as_list(gene_names), r_genes_2_GO, ontology)

def update_go_data(go_data, genes):
    "Update the topGO data with a new set of interesting genes."
    logging.debug('Updating topGO data for %d genes', len(genes))
    return r.updateGoData(go_data, rpy2.robjects.StrVector(as_list(genes)))


def p_value_from_r(p):
    "Caters for results from topGO that are: '< 1e-30'"
    try:
        if p.startswith('<'):
            p = p[1:]
    except AttributeError:
        pass
    return float(p)


def do_go_analysis(go_data, genes, p_value_threshold, method='weight', elim_cut_off=0.01):
    """
    Do a topGO analysis on the genes.

    @arg method: Can be 'classic', 'weight' or 'elim'
    """

    # update the go data with the given genes
    go_data = update_go_data(go_data, genes)

    # get the significant groups
    if 'weight' == method:
        test_stat = r.new(
            "weightCount",
            testStatistic=r.GOFisherTest,
            name="Fisher test",
            sigRatio="ratio"
        )
    elif 'classic' == method:
        test_stat = r.new(
            "classicCount",
            testStatistic=r.GOFisherTest,
            name="Fisher test"
        )
    elif 'elim' == method:
        test_stat = r.new(
            "elimCount",
            testStatistic=r.GOFisherTest,
            name="Fisher test",
            cutOff = elim_cut_off
        )
    else:
        raise ValueError('%s: Unknown topGO method' % method)
    sig_groups = r.getSigGroups(go_data, test_stat)

    args = {method : sig_groups}
    results_unthresholded = r.GenTable(
      go_data,
      #ranksOf="classic",
      orderBy=method,
      **args
    )

    # only keep those results above the threshold
    assert method == results_unthresholded.colnames()[5] # make sure looking at correct column
    # which rows are above the threshold?
    passed = rpy2.robjects.BoolVector([p_value_from_r(p) <= p_value_threshold for p in results_unthresholded[5]])
    return r.subset(results_unthresholded, passed)


def best_p_value(go_analysis):
    "@return: The best p-value in the analysis (assumes sorted)."
    if None != go_analysis and go_analysis.nrow():
        return go_analysis[5][0]
    else:
        return 1.


def yield_stats(go_analysis):
    """
    Extract the information from the go_analysis R object and yield it.

    GO-ID,Term,Annotated,Significant,Expected,weight
    """
    for i in xrange(go_analysis.nrow()):
        yield go_analysis[0][i], go_analysis[1][i], go_analysis[2][i], go_analysis[3][i], p_value_from_r(go_analysis[4][i]), p_value_from_r(go_analysis[5][i])




#set width of console for R printing
logging.info('Setting R console width to %d', options.r_console_width)
set_r_width(options.r_console_width)




if '__main__' == __name__:
    import go

    liver_dev_genes = set((
        "ENSMUSG00000005836",
        "ENSMUSG00000025907",
        "ENSMUSG00000039910",
        "ENSMUSG00000079844",
        "ENSMUSG00000006818",
        "ENSMUSG00000044147",
        "ENSMUSG00000024927",
        "ENSMUSG00000009569",
        "ENSMUSG00000075401",
    ))

    blood_coagulation_genes = set((
        "ENSMUSG00000022149",
        "ENSMUSG00000029664",
        "ENSMUSG00000059481",
        "ENSMUSG00000024909",
        "ENSMUSG00000071311",
        "ENSMUSG00000052229",
        "ENSMUSG00000024386",
        "ENSMUSG00000024940",
        "ENSMUSG00000031443",
    ))

    method = 'weight'
    logging.info('Using topGO method: %s', method)

    all_genes = liver_dev_genes.union(blood_coagulation_genes)

    genes_2_GO = go.get_ensembl_go_annotations(all_genes)

    initialise_topgo()
    r_genes_2_GO = convert_genes_2_GO(genes_2_GO)
    go_data = create_go_data(all_genes, r_genes_2_GO, 'BP')

    p_value_threshold = 5e-2
    set_r_width(240)

    liver_go_analysis = do_go_analysis(go_data, liver_dev_genes, p_value_threshold, method=method)
    logging.info('Analysed liver development genes')
    r_print(liver_go_analysis)

    blood_go_analysis = do_go_analysis(go_data, blood_coagulation_genes, p_value_threshold, method=method)
    logging.info('Analysed blood coagulation genes')
    r_print(blood_go_analysis)

    ps = pretty_str(blood_go_analysis)
    print ps
