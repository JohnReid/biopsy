#
# Copyright John Reid 2006
#

import sys, biopsy

remome_file = sys.argv[1]
analysis_file = sys.argv[2]

print 'Loading from:', remome_file
analysis = biopsy.RemomeAnalysis( biopsy.Remome.load( remome_file ) )

pssms = biopsy.get_transfac_pssm_accessions( biopsy.get_default_transfac_pssm_filter() )
threshold = 0.05
phylo_threshold = 0.02

analysis.analyse_and_serialise(
        analysis_file,
        pssms,
        threshold,
        phylo_threshold )
