#
# Copyright John Reid 2008
#


"""
Code to access INOH data.
"""

from shared import *
import elementtree.ElementTree as ET

bp_ns = "http://www.biopax.org/release/biopax-level2.owl#"
rdf_ns = "http://www.w3.org/1999/02/22-rdf-syntax-ns#"

def xref_value(xref):
    return xref.attrib['{%s}resource' % rdf_ns]

def load_inoh_pathway(filename):
    return ET.parse(filename)

def get_xrefs_in_inoh_pathway(tree):
    return tree.findall('.//{%s}XREF' % bp_ns)

def is_uniprot_ref(value):
    return value.startswith('#UniProt_')

def get_uniprot_id(value):
    return value.split('_')[1]

def get_uniprot_refs_in_inoh_pathway(tree):
    return imap(get_uniprot_id, ifilter(is_uniprot_ref, imap(xref_value, get_xrefs_in_inoh_pathway(tree))))

def get_inoh_pathways():
    pathway_dir = os.path.join(biopsy.get_data_dir(), 'KEGG', 'pathways')
    return [
      (file, set(l.strip() for l in open(os.path.join(pathway_dir, file))))
      for file in os.listdir(pathway_dir)
    ]


def yield_uniprot_gene_mappings(gene_ids):
    import biopsy.identifiers.biomart as biomart, csv
    query = biomart.new_query()
    dataset = biomart.add_dataset(query, 'mmusculus_gene_ensembl')
    biomart.add_filter(dataset, 'ensembl_gene_id', ",".join(gene_ids))
    biomart.add_attribute(dataset, 'ensembl_gene_id')
    biomart.add_attribute(dataset, 'unified_uniprot_accession')
    for row in csv.reader(biomart.execute_query(query), delimiter=','):
        if row[1]:
            yield biopsy.DbRef.parse_as(row[0], biopsy.db.ensembl), biopsy.DbRef.parse_as(row[1], biopsy.db.swissprot)

def get_uniprot_gene_map(gene_ids):
    logging.info('Querying Ensembl biomart for %d genes UNIPROT accessions', len(gene_ids))
    result = cookbook.DictOfSets()
    for ensembl_id, uniprot_id in chain_biomart_queries(yield_uniprot_gene_mappings, gene_ids):
        result[ensembl_id].add(uniprot_id)
    logging.info('Retrieved UNIPROT accessions for %d genes', len(result))
    return result

def reverse_map(m):
    result = cookbook.DictOfSets()
    for key, values in m.iteritems():
        for v in values:
            result[v].add(key)
    return result

if '__main__' == __name__:
    try:
        uniprot_2_ensembl
    except NameError:
        ensembl_2_uniprot = get_uniprot_gene_map(map(str, genes))
        uniprot_2_ensembl = reverse_map(ensembl_2_uniprot)
    uniprot_ids = set(ref.table for ref in uniprot_2_ensembl.keys())

    inoh_dir = os.path.join(biopsy.get_data_dir(), 'INOH')
    pathway_dir = os.path.join(inoh_dir, 'pathways')
    for file in os.listdir(inoh_dir):
        if file.endswith('.owl'):
            pathway_name = os.path.splitext(file)[0]
            tree = load_inoh_pathway(os.path.join(inoh_dir, file))
            pathway_uniprot_ids = set(get_uniprot_refs_in_inoh_pathway(tree))
            intersection = pathway_uniprot_ids.intersection(uniprot_ids)
            genes = set(chain(*[uniprot_2_ensembl[biopsy.DbRef.parse_as(p, biopsy.db.swissprot)] for p in intersection]))
            logging.info('Found %3d elements and %3d genes of %s pathway in our data.', len(intersection), len(genes), file)
            open(os.path.join(pathway_dir, pathway_name), 'w').write('\n'.join(str(g) for g in genes))
