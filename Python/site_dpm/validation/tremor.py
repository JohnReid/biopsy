#
# Copyright John Reid 2009
#

"""
Groups of genes/TFs defined in Hannenhalli TREMOR paper.
"""

cell_cycle = {
        'E2F1' : 'ENSMUSG00000027490',
        'E2F2' : 'ENSMUSG00000007968',
        'E2F3' : 'ENSMUSG00000016477',
        'E2F4' : 'ENSMUSG00000014859',
        'E2F5' : 'ENSMUSG00000027552',
        'E2F6' : 'ENSMUSG00000057469',
        'E2F7' : 'ENSMUSG00000020185',
        'E2F8' : 'ENSMUSG00000046179',
        'CREB' : 'ENSMUSG00000025958',
        'NF-Y' : 'ENSMUSG00000078344',
}

open('TREMOR-cell-cycle-ensembl.txt', 'w').write('\n'.join(cell_cycle.values()))
