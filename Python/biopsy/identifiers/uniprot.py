#
# Copyright John Reid 2007
#

"""
Code to parse UniProt data.

www.ensembl.org/
"""

import gzip, sys, cookbook, biopsy
from . import lazy
T = biopsy.transfac


_uniprot_file = os.path.join(biopsy.get_data_dir(), 'UniProt', 'uniprot_sprot.dat.gz')

def data():
    """
    Returns a handle to the uniprot data.
    """
    return gzip.open(_uniprot_file)

def yield_records(handle):
    """
    Takes a file like handle and separates it into records.
    """
    record = []
    for line in handle:
        if line.startswith('//'):
            if len(record) > 0:
                yield record
                record = []
        else:
            record.append(line.strip())
    if len(record) > 0:
        yield record

_database_names = set()

def parse_reference(ref_str):
    "Parse a reference entry."
    try:
        #print ref_str
        fields = ref_str.split(';')
        db, id = fields[0].strip().lower(), fields[1]
        _database_names.add(db)
        #print db
        if 'ensembl' == db:
            return T.DbRef.parse_as(id, T.db.ensembl)
    except:
        pass; #print "Unexpected error:", sys.exc_info()[1]

def references_for_record(record):
    """
    Takes a record and returns the accessions and the references.
    """
    accessions = set()
    refs = set()
    for line in record:
        if line.startswith('AC  '):
            for acc in line[4:].split(';'):
                acc = acc.strip()
                if acc:
                    accessions.add(T.DbRef.parse_as(acc, T.db.swissprot))
        if line.startswith('DR  '):
            ref = parse_reference(line[4:])
            if ref:
                refs.add(ref)
    return accessions, refs

def references():
    "Yield all (accessions, references) in the UniProt data."
    for record in yield_records(data()):
        yield references_for_record(record)

def accession_ref_map():
    "Build a map from UniProt ids to external references."
    map = cookbook.DictOfSets()
    for accessions, refs in references():
        for acc in accessions:
            for ref in refs:
                map[acc].add(ref)
    return map

_acc2ref_pickle_file = os.path.join(biopsy.get_data_dir(), 'identifiers', 'uniprot', 'acc2ref.pickle')

acc2ref = lazy.PersistableLazyInitialiser(accession_ref_map, _acc2ref_pickle_file)
