#
# Copyright John Reid 2007
#

"""
Code to access and persist Ensembl data.

www.ensembl.org/
"""


import biopsy
T = biopsy.transfac
import biopsy.identifiers.biomart as biomart
import cookbook, cPickle, csv, sys
from . import lazy
from . import mgi


Species = cookbook.NamedTuple(
        "Species",
        "latin short_name abbreviation long_name"
)
def dataset(species):
    "Return the name of the ensembl dataset for the species."
    return '%s_gene_ensembl' % species.long_name
Species.dataset = dataset


species = {
                'AGAP'      : Species('Anopheles_gambiae', '', '', ''),
                'CG'        : Species('Drosophila_melanogaster', '', '', ''),
                'ENSBTAG'   : Species('Bos_taurus', 'cow', '', 'btaurus'),
                'ENSCAFG'   : Species('Canis_familiaris', 'dog', '', 'cfamiliaris'),
                'ENSCING'   : Species('Ciona_intestinalis', '', '', ''),
                'ENSCPOG'   : Species('Cavia_porcellus', '', '', ''),
                'ENSDARG'   : Species('Danio_rerio', 'zebrafish', '', 'drerio'),
                'ENSG'      : Species('Homo_sapiens', 'human', 'hsap', 'hsapiens'),
                'ENSGALG'   : Species('Gallus_gallus', 'chicken', '', 'ggallus'),
                'ENSMODG'   : Species('Monodelphis_domestica', 'opossum', '', 'mdomestica'),
                'ENSMMUG'   : Species('Macaca_mulatta', '', '', ''),
                'ENSMUSG'   : Species('Mus_musculus', 'mouse', 'mmus', 'mmusculus'),
                'ENSORLG'   : Species('Oryzias_latipes', '', '', ''),
                'ENSOCUG'   : Species('Oryctolagus_cuniculus', '', '', ''),
                'ENSPTRG'   : Species('Pan_troglodytes', 'chimp', '', ''),
                'ENSRNOG'   : Species('Rattus_norvegicus', 'rat', 'rnor', 'rnorvegicus'),
                'ENSXETG'   : Species('Xenopus_tropicalis', 'xenopus', '', 'xtropicalis'),
                'GSTENG'    : Species('Tetraodon_nigroviridis', 'tetraodon', '', ''),
                'SINFRUG'   : Species('Takifugu_rubripes', 'fugu', '', 'trubripes'),
}

def get_orthologs(species1, species2):
    "Returns a dict mapping ensembl genes of species 1 to genes of species 2."
    result = cookbook.DictOfSets()

    # build biomart query and execute
    query = biomart.new_query()
    dataset = biomart.add_dataset(query, species[species2].dataset())
    biomart.add_attribute(dataset, 'ensembl_gene_id')
    biomart.add_attribute(dataset, '%s_ensembl_gene' % species[species1].short_name)
    biomart.print_query(query, 'query.xml')
    for row in csv.reader(biomart.execute_query(query), delimiter=','):
        if row[1]:
            ref_1 = biopsy.transfac.DbRef.parse_as(row[0], biopsy.transfac.db.ensembl)
            ref_2 = biopsy.transfac.DbRef.parse_as(row[1], biopsy.transfac.db.ensembl)
            result[ref_2].add(ref_1)

    return result

#
# Orthologs
#
class _Orthologs(dict):
    """
    A dict keyed by a (species1, species2) tuple, e.g. ('ENSG', 'ENSMUSG') that
    maps orthologs from the first species to the second.

    If an entry is missing, goes to biomart to fill it.
    """

    def __missing__(self, key):
        print 'Querying BioMart for Ensembl %s-%s orthologs' % key
        self[key] = get_orthologs(key[0], key[1])
        return self[key]

_orthologs_pickle_file = 'C:/Data/identifiers/ensembl/orthologs.pickle'

orthologs = lazy.PersistableLazyInitialiser(_Orthologs, _orthologs_pickle_file)



#
# Proteins
#
def proteins_for_species(species):
    "Goes to biomart to map all genes in species to Swissprot proteins"
    from . import entrez
    result = cookbook.DictOfSets()

    # build biomart query for swissprot and execute
    query = biomart.new_query()
    dataset = biomart.add_dataset(query, species)
    biomart.add_attribute(dataset, 'ensembl_gene_id')
    biomart.add_attribute(dataset, 'ensembl_transcript_id')
    biomart.add_attribute(dataset, 'uniprot_swissprot_accession')
    for row in csv.reader(biomart.execute_query(query), delimiter=','):
        gene_ref = biopsy.transfac.DbRef.parse_as(row[0], biopsy.transfac.db.ensembl)
        transcript_ref = biopsy.transfac.DbRef.parse_as(row[1], biopsy.transfac.db.ensembl)
        if row[2]:
            protein_ref = biopsy.transfac.DbRef.parse_as(row[2], biopsy.transfac.db.swissprot)
            result[gene_ref].add((transcript_ref, protein_ref))

    # build biomart query for entrez protein and execute
    query = biomart.new_query()
    dataset = biomart.add_dataset(query, species)
    biomart.add_attribute(dataset, 'ensembl_gene_id')
    biomart.add_attribute(dataset, 'ensembl_transcript_id')
    biomart.add_attribute(dataset, 'protein')
    for row in csv.reader(biomart.execute_query(query), delimiter=','):
        gene_ref = biopsy.transfac.DbRef.parse_as(row[0], biopsy.transfac.db.ensembl)
        transcript_ref = biopsy.transfac.DbRef.parse_as(row[1], biopsy.transfac.db.ensembl)
        if row[2]:
            protein_acc = row[2]
            #print "'%s'" % protein_acc
            if protein_acc in entrez.mouse_proteins().acc2id:
                for protein_id in entrez.mouse_proteins().acc2id[protein_acc]:
                    ref = T.DbRef(T.db.entrez_protein, "", protein_id)
                    result[gene_ref].add((transcript_ref, ref))

    return result

class _Proteins(dict):
    """
    A dict keyed by a species, e.g. ('ENSMUSG') that
    maps genes to swissprot proteins.

    If an entry is missing, goes to biomart to fill it.
    """

    species = {
            'ENSMUSG'               : 'mmusculus_gene_ensembl', # mouse
    }

    def __missing__(self, key):
        print 'Querying BioMart for Ensembl mouse proteins'
        species = _Proteins.species[key]
        value = proteins_for_species(species)
        self[key] = value
        return value

_proteins_pickle_file = 'C:/Data/identifiers/ensembl/proteins.pickle'

proteins = lazy.PersistableLazyInitialiser(_Proteins, _proteins_pickle_file)




#
# MGI
#
def _mgi_ids_from_biomart():
    "Goes to biomart to map all mouse genes to MGI ids"
    from . import entrez
    result = cookbook.DictOfSets()

    # build biomart query for swissprot and execute
    query = biomart.new_query()
    dataset = biomart.add_dataset(query, 'mmusculus_gene_ensembl')
    biomart.add_attribute(dataset, 'ensembl_gene_id')
    biomart.add_attribute(dataset, 'external_gene_id')
    for row in csv.reader(biomart.execute_query(query), delimiter=','):
        gene_ref = biopsy.transfac.DbRef.parse_as(row[0], biopsy.transfac.db.ensembl)
        if row[1]:
            mgi_acc = row[1]
            if mgi_acc in mgi.acc2id():
                yield gene_ref, mgi.acc2id()[mgi_acc]

def _mgi_id_map():
    "Returns a dict mapping ensembl refs to mgi ids"
    return dict(_mgi_ids_from_biomart())

_mgi_pickle_file = 'C:/Data/identifiers/ensembl/mgi.pickle'

mgi_ids = lazy.PersistableLazyInitialiser(_mgi_id_map, _mgi_pickle_file)



if '__main__' == __name__:

    mouse_proteins = proteins()['ENSMUSG']
    g = mouse_proteins.keys()[0]
    print mouse_proteins[g]

    human_mouse_orthologs = orthologs()[('ENSG', 'ENSMUSG')]
    print human_mouse_orthologs.keys()[:10]

    persist_all()
