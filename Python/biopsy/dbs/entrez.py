#
# Copyright John Reid 2009
#

"""
Code to interact with Entrez database.
"""


from Bio import Entrez as E

_my_email = 'my.email@biopsy.com'

def get_sequence(chromosome, seq_start, seq_stop, organism='Homo sapiens', build='GRCh37'):
    "@return: Given sequence as fasta string."
    db = 'nucleotide'
    term = '"%s"[Organism] AND "chromosome %s" AND "%s" NOT "contig"' % (organism, chromosome, build)
    search = E.esearch(db=db, term=term, email=_my_email)
    record = E.read(search)
    id = record['IdList'][0]
    record = E.efetch(db=db, id=id, rettype='fasta', seq_start=seq_start, seq_stop=seq_stop)
    return record.read()

if '__main__' == __name__:
    print get_sequence('Y', 2000040, 2000140)
