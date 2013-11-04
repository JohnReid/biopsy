#
# Copyright John Reid 2009
#

"""
Package to handle SymAtlas data stored in the biopsy data directory.
"""

import os, csv, numpy, biopsy, cookbook.cache_decorator
from itertools import imap

datasets = [
    'heart',
    'liver',
]

def average_replicates(replicated_data):
    assert len(replicated_data) % 2 == 0
    return (numpy.array(replicated_data[0::2]) + numpy.array(replicated_data[1::2])) / 2.

def parse_complete_data_set(f):
    "@return: (dataset, tissues, probes, fold_changes)"
    csvreader = csv.reader(f, delimiter='\t')
    headers = csvreader.next()
    #print headers
    expression = list()
    probes = list()
    for row in csvreader:
        assert len(row) == len(headers)
        probes.append(row[0])
        expression.append(average_replicates(map(float, row[1:])))
    return headers[0], headers[1::2], probes, numpy.array(expression)

def get_data_dir():
    return os.path.join(biopsy.get_data_dir(), 'SymAtlas')

def get_data_filename(file):
    return os.path.join(get_data_dir(), file)

def pickled_cached_method(name):
    "Decorator to store output of methods in data directory."
    return cookbook.cache_decorator.pickled_cached_method(os.path.join(get_data_dir(), '%s.pickle' % name))

def annotation_table_filename(dataset):
    return get_data_filename('%s.tsv' % dataset)

def parse_annotation_table(f):
    for row in csv.reader(f, delimiter='\t'):
        (
            name,
            desc,
            taxon,
            reporter_accs,
        ) = row[:4]
        if len(row) > 4:
            (
                genome_locs,
                locuslink_accs,
                refseq_accs,
                unigene_accs,
                uniprot_accs,
                ensembl_accs,
                aliases,
                function_ann,
                protein_family_ann,
                gene_family_ann
            ) = row[4:]
            ensembl_accs = ensembl_accs.split(';')
        else:
            genome_locs = ''
            locuslink_accs = ''
            refseq_accs = ''
            unigene_accs = ''
            uniprot_accs = ''
            ensembl_accs = []
            aliases = ''
            function_ann = ''
            protein_family_ann = ''
            gene_family_ann = ''
        yield (
            name,
            desc,
            taxon,
            reporter_accs,
            genome_locs,
            locuslink_accs,
            refseq_accs,
            unigene_accs,
            uniprot_accs,
            ensembl_accs,
            aliases,
            function_ann,
            protein_family_ann,
            gene_family_ann
        )

def ensembl_genes_in_dataset(dataset):
    import biopsy
    f = open(annotation_table_filename(dataset))
    f.next()
    for row in parse_annotation_table(f):
        for ensembl in row[9]:
            ref = biopsy.DbRef.try_to_parse(ensembl)
            if ref.table == 'ENSMUSG':
                yield ref

def fold_change_above_median(expression, fold_change=10.):
    median = numpy.median(expression, axis=0)
    return expression >= fold_change*median

def fold_change_above_or_below_median(expression, fold_change=10.):
    median = numpy.median(expression, axis=0)
    return (expression >= fold_change*median) + (expression <= median/fold_change)

def match_probe_sets(match_matrix, probes):
    return [
        set(probes[i] for i in numpy.where(tissue_match)[0])
        for tissue_match
        in match_matrix
    ]

def probe_to_ensembl_transcript(
    probe=None,
    chromosome=None,
    refseq=None,
    unigene=None,
    unknown2=None,
    unknown3=None,
    name=None,
    desc=None,
    ensembl_transcript=None,
    unknown4=None,
):
    return probe, ensembl_transcript


@pickled_cached_method('probes-to-transcripts')
def probes_to_transcripts():
    return dict(
        probe_to_ensembl_transcript(*annotation)
        for annotation
        in csv.reader(open(get_data_filename('gnf1m.NEW_ANNO6.tsv')), delimiter='\t')
        if len(annotation) > 8 and annotation[8]
    )

@pickled_cached_method('transcripts-to-genes')
def transcripts_to_genes():
    import biopsy.identifiers.biomart as B
    return dict(
        B.quick_query(
            dataset='mmusculus_gene_ensembl',
            attributes=['ensembl_transcript_id', 'ensembl_gene_id'],
            filters=()
        )
    )

@pickled_cached_method('probes-to-genes')
def probes_to_genes():
    return dict(
        (probe, transcripts_to_genes().get(transcript, None))
        for probe, transcript
        in probes_to_transcripts().iteritems()
    )

@pickled_cached_method('expression-data')
def expression_data():
    "@return: (dataset, tissues, probes, fold_changes)"
    return parse_complete_data_set(open(get_data_filename('GNF1Mdata.txt')))




if '__main__' == __name__:
    heart_genes = set(imap(str, ensembl_genes_in_dataset('heart')))
    print '%d heart genes' % len(heart_genes)
    liver_genes = set(imap(str, ensembl_genes_in_dataset('liver')))
    print '%d liver genes' % len(liver_genes)

    # get the expression data
    dataset, tissues, probes, expression = expression_data()

    # which probes are interesting in which tissues
    fold_change = 10.
    highly_expressed = fold_change_above_or_below_median(expression, fold_change=fold_change)
    highly_expressed_probes = match_probe_sets(highly_expressed.T, probes)
    P.figure()
    P.bar(range(len(tissues)), map(len, highly_expressed_probes))
    P.xlim(max=len(tissues))

    # map the probes onto ensembl genes
    genes = probes_to_genes()
    highly_expressed_genes = [
        set(genes[p] for p in probe_set if p in genes)
        for probe_set
        in highly_expressed_probes
    ]
