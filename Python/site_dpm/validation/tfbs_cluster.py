#
# Copyright John Reid 2009
#

"""
Sets of transcription factors defined by the TFBScluster software.
"""

haematopoietic_tfs = [
    'Jun', # AP1
    'Runx1', # AML1
    'CEBPa', # CEBP
    'EBF1', # EBF
    # 'EBOX',
    'Qrsl1', # EBOX-GATA
    'Myc', # EBOX (c-Myc)
    # 'ETS', # EHF?????
    # 'GATA',
    # 'HMG',
    'Ikzf1', # Ikaros
    'MEF2c', # MEF2
    'MEIS1',
    'MYB',
    # 'NBOX',
    'Nfat5', # NFAT
    # 'NFAT-AP1',
    'Rela', # NFKB
    'Otx1', # OTX
    'PAX5',
    'SP1',
]
"""
Haematopoietic TFBS
"""


liver_tfs = [
    'HNF1a', # HNF1
    'Foxe3', # HNF3
    'HNF4a', # HNF4
    'CEBPa', # CEBP
]
"""
Liver study TFBS (based on work by Krivan and Wasserman, 2001)
"""


muscle_tfs = [
    'MEF2c', # MEF2
    'SP1',
    'SRF',
    'MyoD1', # EBOX (MyoD)
    'TEF',
]
"""
Muscle study TFBS (based on work by Wasserman and Fickett, 1998)
"""

tf_sets = {
    'Haematopoietic' : haematopoietic_tfs,
    'Liver'          : liver_tfs,
    'Muscle'         : muscle_tfs
}

import biopsy.identifiers.synergizer as S
import biopsy.identifiers.biomart as B

for tag, tfs in tf_sets.iteritems():
    translated = S.translate('ensembl', 'Mus musculus', 'mgi_symbol', 'ensembl_gene_id', tfs)
    print translated
    print 'Ensembl translated: %d/%d' % (S.how_many_have_translations(translated), len(tfs))
    ensembl_genes = S.get_translations(translated)
    open('%s-ensembl.txt' % tag, 'w').write('\n'.join(ensembl_genes))
