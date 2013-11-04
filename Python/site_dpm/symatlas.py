#
# Copyright John Reid 2009
#

"""
Code to validate transcriptional programs against SymAtlas data.
"""

from shared import *
import biopsy.data.symatlas as SA
import tp_threshold, gene_set_enrichment

@log_exceptions()
@caching_decorator('symatlas')
def symatlas():
    logging.info('Analysing SymAtlas expression data.')
    probes_to_genes = SA.probes_to_genes()
    dataset, tissues, probes, fold_changes = SA.expression_data()
    highly_expressed = SA.fold_change_above_median(fold_changes, fold_change=options.symatlas_fold_change_threshold)
    highly_expressed_probes = SA.match_probe_sets(highly_expressed.T, probes)
    highly_expressed_genes = [
        set(probes_to_genes[p] for p in hep if p in probes_to_genes and probes_to_genes[p])
        for hep in highly_expressed_probes
    ]
    if False:
        import pylab as P
        P.figure()
        P.bar(range(len(tissues)), map(len, highly_expressed_probes))
        P.xlim(max=len(tissues))
    transcriptional_programs, factor_universe, target_universe = tp_threshold.threshold_tps()
    tester = gene_set_enrichment.TpEnrichmentTester(factor_universe, target_universe, p_value_threshold=1e-3)
    for tissue, tissue_genes in zip(tissues, highly_expressed_genes):
        logging.debug('Testing %d transcriptional programs\' factors for enrichment in genes over-expressed in tissue %s', len(transcriptional_programs), tissue)
        for tp, (test_drawn, test_size, test_complement_size, draws, p_value) in tester.test_transcriptional_program_factors(transcriptional_programs, tissue_genes):
            logging.info(
                'TP:%4d; %4d in program; %4d/%5d over-expressed in % -32s; %4d in intersection; p-value=%e',
                tp.k, draws, test_size, test_complement_size+test_size, tissue, test_drawn, p_value
            )
        for tp, (test_drawn, test_size, test_complement_size, draws, p_value) in tester.test_transcriptional_program_targets(transcriptional_programs, tissue_genes):
            logging.info(
                'TP:%4d; %4d in program; %4d/%5d over-expressed in % -32s; %4d in intersection; p-value=%g',
                tp.k, draws, test_size, test_complement_size+test_size, tissue, test_drawn, p_value
            )

if '__main__' == __name__:
    symatlas()
