#
# Copyright John Reid 2007
#

import biopsy, csv

def transfac_refs( table ):
    "Yield all the references in the transfac table"
    for e in table.all():
        for r in e.db_refs:
            yield r

def transfac_ref_databases( table ):
    "Return the databases linked to from the table"
    dbs = set()
    for r in transfac_refs( table ):
        dbs.add( r.db )
    return dbs

print "\n".join(
        str( db )
        for db
        in transfac_ref_databases( biopsy.transfac.Molecule ) )
raise RuntimeError

def transfac_gene_ensembl_species():
    "Return a set of ensembl species that the transfac gene table references"
    species = set()
    for r in transfac_refs( biopsy.transfac.Gene ):
        if biopsy.db.ensembl == r.db:
            species.add( r.table )
    return species
#print transfac_gene_ensembl_species()


class Mapping( dict ):
    """A mapping from keys to lists of values"""
    def get( self, k ):
        if k in self: return self[k]
        else: return []
    def add( self, k, v ):
        if k in self: self[k].append( v )
        else: self[k] = [ v ]

def csv_reference_mapping( filename, mapped_to_db, mapped_from_db ):
    """Reads a csv file where each line starts with a mapped to reference and
    then has mapped from references"""
    reader = csv.reader( open( filename, "r" ) )
    for row in reader:
        mapped_to = biopsy.DbRef.parse_as( row[0], mapped_to_db )
        for r in row[1:]:
            mapped_from = biopsy.DbRef.parse_as( r, mapped_from_db )
            yield mapped_from, mapped_to

def load_ensembl_mouse_homologs():
    """Loads ensembl mouse homologs into a mapping"""
    mapping = Mapping()
    try:
        for from_ref, to_ref in  csv_reference_mapping(
                "c:/data/identifiers/ensembl_mouse_homologs.csv",
                biopsy.db.ensembl,
                biopsy.db.ensembl
        ):
            mapping.add( from_ref, to_ref )
    except: pass
    return mapping

def count_refs_in_mapping( references, mapping ):
    """Count the number of references in the mapping"""
    count = 0
    for r in references:
        if r in mapping:
            count += 1
    return count

def transfac_references( table ):
    "Yield all the entries in the transfac table as DbRefs"
    for e in table.all():
        # print e.acc.as_db_ref()
        yield e.acc.as_db_ref()


def map_transfac_entry( t, mapping ):
    "Add a transfac entry to the mapping"
    if isinstance( t, biopsy.transfac.Gene ): # a gene
        g_ref = t.acc.as_db_ref()
        for r in t.db_refs:
            if r in mapping:
                mapping.add( g_ref, r )
    elif isinstance( t, biopsy.transfac.Factor ): # a factor
        ref = t.acc.as_db_ref()
        g_ref = t.gene.as_db_ref()
        for e in mapping.get( g_ref ):
            mapping.add( ref, e )
    elif (
            isinstance( t, biopsy.transfac.Matrix )
            or isinstance( t, biopsy.transfac.Site ) # a matrix or site
    ):
        ref = t.acc.as_db_ref()
        for f in t.factors:
            f_ref = f.link.as_db_ref()
            for e in mapping.get( f_ref ):
                mapping.add( ref, e )
    else: raise RuntimeError( "Cannot add %s to mapping" % str( s ) )

def map_transfac_entries( table, mapping ):
    """Maps transfac entries from the given table into the given mapping"""
    for e in table.all():
        map_transfac_entry( e, mapping )
    print "%d %s out of %d mapped" % (
            count_refs_in_mapping( transfac_references( table ), mapping ),
            str( table ),
            len( table.all() )
    ) # how many did we map?

def map_all_transfac_entries( mapping ):
    """Maps all the transfac entries into the given mapping"""
    map_transfac_entries( biopsy.transfac.Gene, mapping )
    map_transfac_entries( biopsy.transfac.Factor, mapping )
    map_transfac_entries( biopsy.transfac.Matrix, mapping )
    map_transfac_entries( biopsy.transfac.Site, mapping )



# load ensembl mouse homologs
mapping = load_ensembl_mouse_homologs()

# add transfac entries to the mapping
map_all_transfac_entries( mapping )
