#
# Copyright John Reid 2008
#

"""
Code to perform a topGO analysis on the transcriptional programs.
"""

import itertools, logging, cookbook, sys
from itertools import imap
from shared import *
import topgo



from rpy2.robjects import r
from rpy2.rinterface import RRuntimeError


def get_ensembl_go_annotations(genes):
    "@return: A map from the given genes to sets of go annotations."
    import biopsy.identifiers.biomart as biomart
    logging.info('Querying Ensembl biomart for GO annotations of %d genes', len(genes))
    result = cookbook.DictOfLists()
    for id_attr, evidence_attr in [
      ('go_biological_process_id', 'go_biological_process_linkage_type'),
      ('go_cellular_component_id', 'go_cellular_component_linkage_type'),
      ('go_molecular_function_id', 'go_molecular_function_linkage_type'),
    ]:
        query = biomart.new_query()
        dataset = biomart.add_dataset(query, 'mmusculus_gene_ensembl')
        biomart.add_attribute(dataset, 'ensembl_gene_id')
        biomart.add_attribute(dataset, id_attr)
        biomart.add_attribute(dataset, evidence_attr)
        filter = biomart.add_filter(dataset, name='ensembl_gene_id', value='')
        for chunk in biomart.split_big_list((str(g) for g in genes), 50):
            #logging.info('Querying Ensembl biomart for chunk of %d genes', len(chunk))
            filter.set('value', ','.join(chunk))
            for row in biomart.yield_csv_query_results(query):
                if row[2] not in options.go_evidence_codes_to_ignore:
                    result[row[0]].append(row[1])
    logging.info('Found %d go annotations', sum(len(v) for v in result.values()))
    return result


@global_cached_method('all-go-annotations')
def get_all_ensembl_go_annotations():
    "@return: A map from ensembl genes to sets of go annotations."
    import biopsy.identifiers.biomart as biomart
    logging.info('Querying Ensembl biomart for all GO annotations')
    result = cookbook.DictOfLists()
    for id_attr, evidence_attr in [
        ('go_biological_process_id', 'go_biological_process_linkage_type'),
        ('go_cellular_component_id', 'go_cellular_component_linkage_type'),
        ('go_molecular_function_id', 'go_molecular_function_linkage_type'),
    ]:
        for row in biomart.quick_query(
            dataset='mmusculus_gene_ensembl',
            attributes=['ensembl_gene_id', id_attr, evidence_attr]
        ):
            if row[2] not in options.go_evidence_codes_to_ignore and row[1]:
                result[row[0]].append(row[1])
    logging.info('Found %d go annotations', sum(len(v) for v in result.values()))
    return result


def prepare_topGO(factor_universe, target_universe):
    #
    # Load interface to topGO R package
    #
    import biopsy
    topgo.initialise_topgo()
    def db_ref_as_r(ref):
        return str(ref)
    biopsy.DbRef.as_r = db_ref_as_r
    genes_2_GO = get_all_ensembl_go_annotations()
    logging.info(
        '%d / %d factors have no mapping to GO',
        len([f for f in factor_universe if f not in genes_2_GO]),
        len(factor_universe)
    )
    for f in factor_universe: # ensure they have empty entries
        genes_2_GO[f]
    logging.info(
        '%d / %d targets have no mapping to GO',
        len([t for t in target_universe if t not in genes_2_GO]),
        len(target_universe)
    )
    for t in target_universe: # ensure they have empty entries
        genes_2_GO[t]
    r_genes_2_GO = topgo.convert_genes_2_GO(genes_2_GO)
    return r_genes_2_GO






def p_value_from_r(p):
    if p.startswith('<'):
        p = p[1:]
    return float(p)

def test_go_enrichment(genes, go_data, p_value_threshold, ontology, output_dir):
    """
    Test if go enrichment produces significant p-values for random sets of genes
    """
    from random import sample
    return r.goEnrichment(
      go_data,
      sample(genes, 70),
      p_value_threshold,
      os.path.join(output_dir, 'enrichment-test')
    )
# test_go_enrichment(factors, factors_go_data, p_value_threshold, ontology, output_dir); raise RuntimeError()




def threshold_results(results, threshold):
    if not results.classic:
        return results
    return r.subset(results, [c < threshold for c in imap(p_value_from_r, results.classic)])




class GoContext(object):
    """
    Global data/context for GO analysis.
    """

    def __init__(self, factor_universe, target_universe, factors_go_data, targets_go_data):
        self.factor_universe = factor_universe
        self.target_universe = target_universe
        self.factors_go_data = factors_go_data
        self.targets_go_data = targets_go_data




def initialise_go_context(factor_universe, target_universe, go_ontologies):
    """
    Global initialisation of GO context/data.
    """
    genes_2_GO = prepare_topGO(factor_universe, target_universe)
    factors_go_data = dict(
        (ontology, topgo.create_go_data(factor_universe, genes_2_GO, ontology))
        for ontology in go_ontologies
    )
    targets_go_data = dict(
        (ontology, topgo.create_go_data(target_universe, genes_2_GO, ontology))
        for ontology in go_ontologies
    )
    go_context = GoContext(factor_universe, target_universe, factors_go_data, targets_go_data)
    return genes_2_GO, go_context



def try_go_analysis(go_data, interesting, p_value_threshold, topgo_method):
    try:
        return topgo.do_go_analysis(go_data, interesting, p_value_threshold, method=topgo_method)
    except RRuntimeError:
        print sys.exc_info()
        logging.warning('Could not run GO analysis')
        return None




class TPGoAnalysis(object):
    """
    The GO analysis of a particular transcriptional program.
    """

    def __init__(self, tp, go_context, p_value_threshold, topgo_method):
        logging.info('Calculating GO enrichment analysis for transcriptional program; %d' % tp.k)
        self.tp = tp
        self.go_context = go_context
        self.p_value_threshold = p_value_threshold
        self.topgo_method = topgo_method
        self.factors_go_analysis = dict(
            (
                ontology,
                try_go_analysis(
                    go_data,
                    tp.factors,
                    p_value_threshold,
                    topgo_method
                )
            )
            for ontology, go_data
            in self.go_context.factors_go_data.iteritems()
        )
        self.targets_go_analysis = dict(
            (
                ontology,
                try_go_analysis(
                    go_data,
                    tp.targets,
                    p_value_threshold,
                    topgo_method
                )
            )
            for ontology, go_data
            in self.go_context.targets_go_data.iteritems()
        )


    def print_go_analysis(self, f, tag, number, ontology, threshold, go_analysis, log=False):
        if None != go_analysis and go_analysis.nrow():
            output = '************ TP: %d **************; # %s=%d; ontology=%s\n%s' % (self.tp.k, tag, number, ontology, topgo.pretty_str(go_analysis))
            print >> f, output
            if log:
                logging.info(output)


    def print_go_analyses(self, f, threshold, log=False):
        for ontology, go_analysis in self.factors_go_analysis.iteritems():
            self.print_go_analysis(f, 'factors', len(self.tp.factors), ontology, threshold, go_analysis, log)
        for ontology, go_analysis in self.targets_go_analysis.iteritems():
            self.print_go_analysis(f, 'targets', len(self.tp.targets), ontology, threshold, go_analysis, log)
