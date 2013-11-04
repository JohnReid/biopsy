#
# Copyright John Reid 2010
#

"""
Locate nearest gene for a genomic region.
"""

import os, itertools, logging, gzip, csv
import biopsy.data, cookbook.dicts as D 


gff_files = {
    'Drosophila_melanogaster' : {
        'dmel_r4.1' : 'dmel-%s-r4.1.1.gff.gz',
    },
}

cpg_files = {
    'Drosophila_melanogaster' : {
        'dmel_r4.1' : 'dmel-%s-chromosome-r4.1.1.cpgreport',
    },
}

chromosomes = {
    'Drosophila_melanogaster' : ('2L', '2R', '3L', '3R', '4', 'X'),
}

def data_dir():
    "@return: The data directory."
    return os.path.join(biopsy.data.data_dir(), 'FlyBase')

def genes_dir():
    "@return: The genes directory."
    return os.path.join(data_dir(), 'genes')

def genomes_dir():
    "@return: The genomes directory."
    return os.path.join(data_dir(), 'genomes')

def fasta_dir(organism, release):
    "@return: The directory of fasta files."
    return os.path.join(genomes_dir(), organism, release, 'fasta')

def cpg_dir(organism, release):
    "@return: The directory of CpG files."
    return os.path.join(genomes_dir(), organism, release, 'fasta', 'CpG')

def cpg_file(organism, release, chr):
    "@return: The CpG filename."
    return os.path.join(cpg_dir(organism, release), cpg_files[organism][release] % chr)

def gff_dir(organism, release):
    "@return: The directory of GFF files."
    return os.path.join(genomes_dir(), organism, release, 'gff')

def gff_file(organism, release, chr):
    "@return: The GFF filename."
    return os.path.join(gff_dir(organism, release), gff_files[organism][release] % chr)

def gene_map_file():
    "@return: The name of the file mapping genes to locations."
    return os.path.join(genes_dir(), 'gene_map_table_fb_2010_03.tsv')

def comment(line):
    "@return: True if the line is a comment."
    return line.startswith('#')

def parse_seq_loc(loc):
    "Parse '3R:5513119..5518516(-1)' into id, start, end, orientation"
    id_, pos = loc.split(':')
    pos, orientation = pos.split('(')
    orientation = orientation[:-1]
    start, end = map(int, pos.split('..'))
    return id_, start, end, orientation

def parse_gene_map_file(f):
    "@return: Yield entries in gene map file."
    for row in csv.reader(itertools.ifilterfalse(comment, f), delimiter='\t'):
        if row:
            symbol, primary_FBid, recomb_loc, cytogenetic_loc, seq_loc = row
            if seq_loc:
                seq_loc = parse_seq_loc(seq_loc)
            yield symbol, primary_FBid, recomb_loc, cytogenetic_loc, seq_loc

def build_melanogaster_gene_map():
    "Build an interval map for each melanogaster chromosome."
    from pyitl import IntIntervalMap, IntInterval
    genes = D.DictOf(IntIntervalMap)
    for symbol, primary_FBid, recomb_loc, cytogenetic_loc, seq_loc in parse_gene_map_file(open(gene_map_file(), 'rb')):
        if -1 == symbol.find('\\'): # assume anything without \\ is melanogaster. Is this right?
            if seq_loc:
                id_, start, end, orientation = seq_loc
                genes[id_].add(IntInterval(start, end), ((symbol, primary_FBid),))
    return genes

def load_gff_file(organism, release, chr, types=None):
    "@return: The features in the GFF file."
    filename = gff_file(organism, release, chr)
    logging.debug('Loading GFF file: %s', filename)
    for record in csv.reader(gzip.open(filename, 'rb'), delimiter='\t'):
        if len(record) and not record[0].startswith('#'):
            if None == types or record[2] in types:
                yield record

    
def load_gff_files(organism, release, *load_gff_args, **load_gff_kwargs):
    "@return: A dict mapping chromosomes to gff features."
    return dict(
        (chr, list(load_gff_file(organism, release, chr, *load_gff_args, **load_gff_kwargs)))
        for chr in chromosomes[organism] 
    )

def interval_maps_from_features(features):
    "@return: An interval map built from the features."
    from pyitl import IntIntervalSet, IntInterval
    result = D.DictOf(IntIntervalSet)
    for f in features:
        result[f[2]].add(IntInterval(int(f[3]), int(f[4])))
    return result

def load_CpG_islands(organism, release, chr):
    "@return: Yield the CpG islands."
    from Bio import SeqIO
    filename = cpg_file(organism, release, chr)
    logging.debug('Loading CpG islands from %s', filename)
    for l in open(filename, 'r'):
        if l.startswith('FT   CpG island'):
            location = l.split()[3]
            yield map(int, location.split('..'))

def interval_map_from_locations(locations):
    "@return: An interval map of the (start, end) tuples in locations."
    from pyitl import IntIntervalSet, IntInterval
    result = IntIntervalSet()
    for start, end in locations:
        result.add(IntInterval(start, end))
    return result
    

if '__main__' == __name__:
    logging.basicConfig(level=logging.DEBUG)
    
    #exons = load_gff_files('Drosophila_melanogaster', 'dmel_r4.1', types=set('exon'))
    #feature_maps = interval_maps_from_features(load_gff_file('Drosophila_melanogaster', 'dmel_r4.1', '2L'))
    
    CpG_islands = interval_map_from_locations(load_CpG_islands('Drosophila_melanogaster', 'dmel_r4.1', '2L'))
    for interval in CpG_islands:
        logging.info(interval)
    
#    genes = build_melanogaster_gene_map()
#    for id_, intervals in genes.iteritems():
#        print id_, intervals.iterative_size()
