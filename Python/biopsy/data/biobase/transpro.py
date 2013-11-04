#
# Copyright John Reid 2009
#

"""
Code to parse biobase TRANSFAC TransPro data files.
"""


from itertools import imap
import os, biopsy, cookbook.cache_decorator

def pickled_cached_method(name):
    "Decorator to store output of methods in data directory."
    return cookbook.cache_decorator.pickled_cached_method(
      os.path.join(get_data_dir(),
      '%s.pickle' % name)
    )


def get_data_dir():
    return os.path.join(biopsy.get_data_dir(), 'biobase', 'transpro')

def transpro_filename():
    return "transpro.dat"

def transpro_file():
    return os.path.join(get_data_dir(), transpro_filename())

def split_embl_line(line):
    return line[:2], line[4:]

class Record(object):
    def __init__(self):
        self.sequence = ''
        self.refs = []

dispatcher = {}

def handle_accession(record, data):
    record.acc = data
dispatcher['AC'] = handle_accession

def handle_identifier(record, data):
    record.id = data
dispatcher['ID'] = handle_identifier

def handle_description(record, data):
    record.description = data
dispatcher['DE'] = handle_description

def handle_gene(record, data):
    record.gene = data
dispatcher['GS'] = handle_gene

def handle_synonyms(record, data):
    record.synonyms = map(str.strip, data.split(','))
dispatcher['SY'] = handle_synonyms

def handle_species(record, data):
    record.species = data
dispatcher['OS'] = handle_species

def handle_database_ref(record, data):
    record.refs.append(data)
dispatcher['DR'] = handle_database_ref

def handle_sequence(record, data):
    record.sequence += data
dispatcher['SQ'] = handle_sequence

def handle_chromosome(record, data):
    record.chromosome = data
dispatcher['CH'] = handle_chromosome

def handle_chromosomal_location(record, data):
    if data.startswith('Sequence Fragment'):
        record.location = data
    elif data.startswith('Score Points'):
        record.score = data
    else:
        raise RuntimeError('Could not parse chromosomal location: %s' % data)
dispatcher['CC'] = handle_chromosomal_location

def build_records(input):
    record = None
    for type_code, data in input:
        # print type_code, data
        if '//' == type_code:
            if record:
                yield record
            record = Record()
        elif type_code in dispatcher:
            dispatcher[type_code](record, data)
    if record:
        yield record

def yield_transpro_records(f):
    for record in build_records(imap(split_embl_line, imap(str.strip, f))):
        yield record

@pickled_cached_method('mouse-promoters')
def get_mouse_promoters():
    return [
        p
        for p
        in yield_transpro_records(open(transpro_file()))
        if hasattr(p, 'species') and p.species == 'mouse, Mus musculus'
    ]

if '__main__' == __name__:
    mouse_promoters = get_mouse_promoters()
