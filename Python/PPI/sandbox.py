#
# Copyright John Reid 2007
#


import biopsy
import biopsy.psimi2 as psimi
psimi = reload(biopsy.psimi2)
import biopsy.identifiers.build_map as build_map
import biopsy.identifiers.entrez as entrez
import biopsy.identifiers.ensembl as ensembl
from cookbook import DictOfSets
T = biopsy.transfac


def not_found(found):
    return [m for m in T.Matrix.all() if m.name.startswith('V') and m.acc.as_db_ref() not in found]

def names_for_factor(f):
    names = set()
    names.add(f.name)
    for s in f.synonyms:
        names.add(s)
    if f.gene:
        names.add(f.gene.entry.name)
    for su in f.subunits:
        names.update(names_for_factor(su.entry))
    return names

def names_for_matrix(m):
    names = set()
    for fl in m.factors:
        f = fl.link.entry
        names.update(names_for_factor(f))
    return names
#
# Test how many matrices get linked to protein interactions
#
if True:
    try:
        matrix_links
        matrix_names
    except NameError:
        print 'Getting database references and names for all matrices'
        map = build_map.build_map()
        matrix_links = dict((m.acc.as_db_ref(), map.links(m.acc.as_db_ref())) for m in T.Matrix.all())
        matrix_names = dict((m.acc.as_db_ref(), names_for_matrix(m)) for m in T.Matrix.all())

    # look at one matrix in particular
    mat = T.Matrix(1151)
    mat_links = map.links(mat.acc.as_db_ref())
    myod_gene = T.DbRef.parse_as('17927', T.db.entrez_gene)
    bind = psimi.networks['BIND']
    myod_in_bind = bind.nodes_for_ref(myod_gene)
    myod_node = bind.nodes_for_name('MYOD')[0]
    myod_interactor = bind.interactor_for_node(myod_node)

    # those pssms we are interested in
    pssm_accs = [
            T.TableLink(acc)
            for acc
            in biopsy.get_transfac_pssm_accessions(
                    biopsy.get_default_transfac_pssm_filter()
            )
    ]

    for db_name in psimi.dbs():
        print
        print db_name
        network = psimi.networks[db_name]
        print network.summary()
        transfac_2_network = psimi.Transfac2Network(map, network)
        print transfac_2_network.coverage_summary(pssm_accs)
