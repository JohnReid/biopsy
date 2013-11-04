#
# Copyright John Reid 2009
#

"""
Code that deals with sets of genes given in the literature.
"""

import transfac_factors, logging
from analysis import get_pssm_to_ensembl_map_min_range
from itertools import imap, chain


liver_ensembl_targets = {
    'G6PC'    : 'ENSMUSG00000078650',
    'IGF1'    : 'ENSMUSG00000020053',
    'PAH'     : 'ENSMUSG00000020051',
    'IGFBP1'  : 'ENSMUSG00000020429',
    'CFB'     : 'ENSMUSG00000024371',
    'FABP2'   : 'ENSMUSG00000023057',
    'GUCA2B'  : 'ENSMUSG00000032978',
    'HOXA4'   : 'ENSMUSG00000000942',
    'SLC34A1' : 'ENSMUSG00000021490',
}
"Targets of liver specific transcription factors in Wasserman's predictive model for liver paper"


pastaa_matrices = {
    'muscle'    : [ 'SRF_Q5_01', 'SRF_01', 'SRF_Q5_02', 'SRF_C', 'MTATA_B' ],
    'heart'     : [ 'MEF2_Q6_01', 'SRF_C', 'RSRFC4_01', 'MTATA_B', 'MEF2_02' ],
    'liver'     : [ 'HNF4_Q6_01', 'HNF1_01', 'HNF4_01', 'HNF1_Q6', 'HNF1_C', ],
    'retina'    : [ 'CRX_Q4', 'CHX10_01' ],
    'leukocyte' : [ 'NFKAPPAB65_01', 'NFKAPPAB_01', 'NFKB_Q6_01', 'CREL_01', 'ETS_Q6' ]
}
"""
Sets of TFs defined in PASTAA paper (Vingron)
"""





tf_sets_by_name = {

    # The transcription factors in "Genome-Wide In Silico Identification of Transcriptional
    # Regulators Controlling the Cell Cycle in Human Cells" Elkon et al. Table 1
    'elkon-cell-cycle-small' : [
        'E2F',
        'NF-Y',
        'NRF-1',
        'CREB',
    ],

    # The transcription factors in "Genome-Wide In Silico Identification of Transcriptional
    # Regulators Controlling the Cell Cycle in Human Cells" Elkon et al. Table 3
    'elkon-cell-cycle' : [
        'E2F',
        'Sp1',
        'NF-Y',
        'NRF-1',
        'ETF',
        'ATF',
        'CREB',
        'Arnt',
        'YY1',
    ],

    #"Liver specific transcription factors as listed in TREMOR paper."
    'tremor-liver' : [
        'HNF1',
        'HNF3',
        'HNF4',
        'cEBP',
    ],

    # Factors interacting according to Luscombe paper
    'luscombe-1' : [
        'SRF',
        'NFAT',
        'HERP1',
    ],

    # Factors interacting according to Luscombe paper
    'luscombe-2' : [
        'SRF',
        'NKX3-1'
    ],

    # Muscle specific transcription factors as listed in TREMOR paper.
    'tremor-muscle' : [
        'MYOD',
        'SRF',
        'MEF2',
        'TEF1',
        'SP1',
    ],

    # Cell-cycle specific transcription factors as listed in TREMOR paper.
    'tremor-cell-cycle' : [
        'E2F1',
        'E2F2',
        'E2F3',
        'E2F4',
        'E2F5',
        'E2F6',
        'E2F7',
        'E2F8',
        'CREB',
        'NF-Y',
    ],

    # Factors involved in muscle-specific transcription according to PASTAA paper
    # Roider-PASTAA-identifying transcription factors associated with sets of co-regulated genes
    #'pastaa-muscle' : [
    #    'SRF',
    #    'MTATA'
    #],

    # Factors involved in heart-specific transcription according to PASTAA paper
    # Roider-PASTAA-identifying transcription factors associated with sets of co-regulated genes
    'pastaa-heart' : [
        'SRF',
        'MEF2',
        'MTATA'
    ],

    # Factors involved in liver-specific transcription according to PASTAA paper
    # Roider-PASTAA-identifying transcription factors associated with sets of co-regulated genes
    'pastaa-liver' : [
        'HNF6',
        'HNF4',
        'HNF1'
    ],

    # Factors involved in retina-specific transcription according to PASTAA paper
    # Roider-PASTAA-identifying transcription factors associated with sets of co-regulated genes
    'pastaa-retina' : [
        'CRX',
        'CHX10'
    ],

    # Factors involved in leukocyte-specific transcription according to PASTAA paper
    # Roider-PASTAA-identifying transcription factors associated with sets of co-regulated genes
    'pastaa-leukocyte' : [
        'NFKB1',
        'c-rel',
        'Elf-1'
    ],

    # Factors involved in haematopoietic-specific transcription according to TFBScluster software
    'tfbs-cluster-haematopoietic' : [
        'Jun', # AP1
        'Runx1', # AML1
        'CEBPa', # CEBP
        'EBF1', # EBF
        'EBOX',
        'Qrsl1', # EBOX-GATA
        'Myc', # EBOX (c-Myc)
        'ETS', # EHF?????
        'GATA',
        'HMG',
        'Ikzf1', # Ikaros
        'MEF2c', # MEF2
        'MEIS1',
        'MYB',
        'NBOX',
        'Nfat5', # NFAT
        'NFAT-AP1',
        'Rela', # NFKB
        'Otx1', # OTX
        'PAX5',
        'SP1',
    ],

    # Factors involved in liver-specific transcription according to TFBScluster software
    'tfbs-cluster-liver' : [
        'HNF1a', # HNF1
        'Foxe3', # HNF3
        'HNF4a', # HNF4
        'CEBPa', # CEBP
    ],

    # Factors involved in muscle-specific transcription according to TFBScluster software
    'tfbs-cluster-muscle' : [
        'MEF2c', # MEF2
        'SP1',
        'SRF',
        'MyoD1', # EBOX (MyoD)
        'TEF',
    ],
}


def map_factors_to_ensembl(name, factor_names, pssm_map):
    result = dict(
        (factor_name, transfac_factors.ens_genes_for(factor_name, pssm_map))
        for factor_name in factor_names
    )
    logging.info('%s\n%s\n', name, '\n'.join('% 20s: %s' % (factor_name, ','.join(imap(str, ensembl))) for factor_name, ensembl in result.iteritems()))
    return result


def extract_ensembl(factor_mapping):
    return set(chain(*[ensembl for ensembl in factor_mapping.values()]))


def map_factor_sets(factor_sets):
    pssm_map = get_pssm_to_ensembl_map_min_range()
    return dict(
        (name, extract_ensembl(map_factors_to_ensembl(name, factor_names, pssm_map)))
        for name, factor_names
        in factor_sets.iteritems()
    )


if '__main__' == __name__:
    tremor_liver_tfs = tf_sets_by_name['tremor-liver']
    pssm_map = get_pssm_to_ensembl_map_min_range()
    tremor_liver_ensembl = map_factors_to_ensembl('tremor-liver', tremor_liver_tfs, pssm_map)
