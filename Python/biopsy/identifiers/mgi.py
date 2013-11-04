#
# Copyright John Reid 2007
#

"""
Code to do with MGI identifiers.
"""


import biopsy
T = biopsy.transfac
import csv, os
from . import lazy

def ensembl_mouse_genes_in_transfac():
    "Returns all the ensembl mouse gene references in transfac gene table."
    mouse_genes = set()
    for g in T.Gene.all():
        for r in g.db_refs:
            if T.db.ensembl == r.db and r.table == 'ENSMUSG':
                mouse_genes.add(r)
    return mouse_genes

def write_ensembl_mouse_genes_file(filename):
    """
    Writes all the ensembl mouse gene references in transfac gene table to a
    file for use at
    http://www.informatics.jax.org/javawi2/servlet/WIFetch?page=batchQF.
    """
    f = open(filename, "w")
    for g in ensembl_mouse_genes_in_transfac():
        f.write(str(g))
        f.write('\n')
    f.close()

accession_map_filename = os.path.join(biopsy.get_data_dir(), 'identifiers', 'mgi', 'MRK_Dump1.rpt')

def parse_accession_map():
    """
    Parse the MGI accession map flat file and yield tuples

    (mgi identifier ref, mgi accession)
    """
    for id, acc in csv.reader(open(accession_map_filename, 'rb'), delimiter='\t'):
        yield (T.DbRef.parse_as(id, T.db.mgi), acc)

def build_acc_2_id_map():
    """
    Returns a dict mapping MGI accessions to ids.
    """
    return dict((acc, id) for id, acc in parse_accession_map())

accession_map_pickle_file = os.path.join(biopsy.get_data_dir(), 'identifiers', 'mgi', 'acc2id.pickle')

acc2id = lazy.PersistableLazyInitialiser(build_acc_2_id_map, accession_map_pickle_file)



mgi_data_dir = os.path.join(biopsy.get_data_dir(), 'MGI')
mrk_ensembl_filename = os.path.join(mgi_data_dir, 'MRK_ENSEMBL.rpt')
gene_association_filename = os.path.join(mgi_data_dir, 'gene_association.mgi')

def build_mgi_ensembl_map():
    """
    @return: A bidirectional-map between MGI identifiers and Ensembl identifiers.
    """
    import csv, cookbook
    from biopsy import DbRef, db
    reader = csv.reader(open(mrk_ensembl_filename, "r"), delimiter='\t')
    result = cookbook.BidirectionalMap()
    for row in reader:
        result.add(DbRef.parse_as(row[0], db.mgi), DbRef.parse_as(row[5], db.ensembl))
    return result

def build_mgi_name_map():
    """
    @return: A map between MGI identifiers and marker names.
    """
    import csv, cookbook
    from biopsy import DbRef, db
    reader = csv.reader(open(mrk_ensembl_filename, "r"), delimiter='\t')
    return dict(
            (DbRef.parse_as(row[0], db.mgi), row[1])
            for row
            in reader
    )

def build_mgi_go_map():
    """
    @return: A map between MGI identifiers and GO ontologies.
    """
    import csv, cookbook
    from biopsy import DbRef, db
    reader = csv.reader(open(gene_association_filename, "r"), delimiter='\t')
    result = cookbook.DictOfSets()
    for row in reader:
        if row[0].startswith('!'):
            continue
        result[DbRef.parse_as(row[1], db.mgi)].add(row[4])
    return result
