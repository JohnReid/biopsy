#
# Copyright John Reid 2007
#


"""
Code to access NCBI Entrez databases
"""

import biopsy
T = biopsy.transfac
import cookbook, cPickle, time
from . import lazy
import csv, cookbook, gzip

def protein_ref(id):
    return T.DbRef(T.db.entrez_protein, "", int(id))

def gene_ref(id):
    return T.DbRef(T.db.entrez_gene, "", int(id))

def parse_ref(ref_text):
    try:
        return T.DbRef.parse(ref_text)
    except:
        try:
            db_name, acc = ref_text.split(':')
            if 'CCDS' == db_name:
                return
            elif 'GOA' == db_name:
                return
            elif 'PDB' == db_name:
                return
            elif 'InterPro' == db_name:
                return
            elif 'GeneID' == db_name:
                return T.DbRef.parse_as(acc, T.db.entrez_gene)
            elif 'UniProtKB/Swiss-Prot' == db_name:
                return T.DbRef.parse_as(acc, T.db.swissprot)
        except:
            print 'Could not parse %s' % ref_text

def mouse_proteins_result():
    from Bio.EUtils import HistoryClient
    search_term = 'mouse[orgn]'
    #search_term = 'MYOD1[Gene name] AND mouse[orgn]' # for testing

    # get a handle to the results
    client = HistoryClient.HistoryClient()
    return client.search(db='protein', term=search_term)

def refs_from_protein_xml(xml):
    import ElementTree.ElementTree as ET
    tree = ET.parse(xml)
    root = tree.getroot()
    for seq in root.findall('GBSeq'):
        seq_acc = seq.find('GBSeq_primary-accession').text
        for other_id in seq.findall('GBSeq_other-seqids/GBSeqid'):
            s = other_id.text.split('|')
            if 'gi' == s[0]:
                seq_id = s[1]
                break
        else:
            continue
        refs = set()
        for feature in seq.findall('GBSeq_feature-table/GBFeature'):
            key = feature.find('GBFeature_key')
            if None != key and 'CDS' == key.text:
                for qualifier in feature.findall('GBFeature_quals/GBQualifier'):
                    name = qualifier.find('GBQualifier_name')
                    if None != name and 'db_xref' == name.text:
                        ref = parse_ref(qualifier.find('GBQualifier_value').text)
                        if ref:
                            refs.add(ref)
        yield seq_acc, T.DbRef.parse_as(seq_id, T.db.entrez_protein), refs

def refs_for_mouse_protein_accs(chunk_size=500):
    "Yield the references for each mouse protein accession in entrez"
    results = mouse_proteins_result()
    results_size = len(results)
    print '# mouse proteins from entrez: %d' % results_size

    for start in xrange(0, results_size, chunk_size):
        size = min(chunk_size, results_size-start)
        end = start+size
        results.retstart = start
        results.retmax = end
        print 'Getting %d->%d' % (start, end)
        for acc, id, refs in refs_from_protein_xml(results.efetch()):
            yield acc, id, refs

ProteinMap = cookbook.NamedTuple('ProteinMap', 'acc_2_id xrefs')

def get_protein_map():
    """
    Query Entrez to get a map from its protein accessions to ids and xrefs
    """
    result = ProteinMap(
            acc_2_id = cookbook.DictOfSets(),
            xrefs = cookbook.DictOfSets()
    )
    for acc, id, refs in refs_for_mouse_protein_accs():
        result.acc_2_id[acc].add(id.acc)
        for ref in refs:
            result.xrefs[id.acc].add(ref)
    return result




_proteins_pickle_file = os.path.join(biopsy.get_data_dir(), 'identifiers', 'entrez', 'proteins.pickle')

proteins = lazy.PersistableLazyInitialiser(get_protein_map, _proteins_pickle_file)

def write_mouse_protein_ids(filename):
    from Bio.EUtils import HistoryClient
    f = open(filename, 'w')
    results = HistoryClient.HistoryClient().search(db='protein', term='mouse[orgn]')
    for id in results.dbids.ids:
        f.write(id)
        f.write('\n')
    f.close()




if False: # no longer used
    def get_protein_accession_map():
        """
        Query Entrez to get a map from its protein accession to ids
        """
        from Bio.EUtils import HistoryClient
        search_term = 'mouse[orgn]'
        #search_term = 'MYOD[Gene name] AND mouse[orgn]'

        # get a handle to the results
        client = HistoryClient.HistoryClient()
        results = client.search(db='protein', term=search_term)
        results_size = len(results)
        dbids = results.dbids.ids
        print '# results: %d' % results_size

        # download them bit by bit
        acc_2_id = cookbook.DictOfSets()
        step = 10000
        for start in xrange(0, results_size, step):
            size = min(step, results_size-start)
            end = start+size
            results.retstart = start
            results.retmax = end
            print 'Getting %d->%d' % (start, end)
            for id, acc in zip(
                    dbids[start:end],
                    results.efetch(retmode='text',rettype='acc')
            ):
                acc = acc.strip().split('.')[0]
                print acc, id
                acc_2_id[acc].add(
                        biopsy.transfac.DbRef.parse_as(id, biopsy.transfac.db.entrez_protein))

        return acc_2_id

    def get_protein_accession_map_2():
        import elementtree.ElementTree as ET
        search_term = 'mouse[orgn]'
        #search_term = 'MYOD[Gene name] AND mouse[orgn]'

        client = HistoryClient.HistoryClient()
        results = client.search(db='protein', term=search_term)
        results_size = len(results)
        print '# results: %d' % results_size

        acc_2_ids = cookbook.DictOfSets()
        step = 5000
        for start in xrange(0, results_size, step):
            results.retstart = start
            results.retmax = min(step, results_size-start)
            #results.retmax = 1000
            start = time.time()
            summary = results.summary()
            print 'Retrieving summary: %f secs' % (time.time() - start)
            for entry in summary:
                acc = entry.dataitems['Caption'].encode().strip().split('.')[0]
                acc_2_ids[acc].add(
                        biopsy.transfac.DbRef.parse_as(
                                entry.id.encode(),
                                biopsy.transfac.db.entrez_protein))




Gene = cookbook.NamedTuple('Gene', 'id symbol refs')

def parse_gene_info(f, taxid_filter='10090'):
    """
    Parse a gene info file object.

    E.g. ftp://ftp.ncbi.nih.gov/gene/DATA/gene_info.gz

    The format as defined by NCBI:
    ===========================================================================
    gene_info
    ---------------------------------------------------------------------------
               tab-delimited
               one line per GeneID
               Column header line is the first line in the file.
    ---------------------------------------------------------------------------

    tax_id:
               the unique identifier provided by NCBI Taxonomy
               for the species or strain/isolate

    GeneID:
               the unique identifier for a gene
               ASN1:  geneid
               --note:  for genomes previously available from LocusLink,
                        the identifiers are equivalent

    Symbol:
               the default symbol for the gene
               ASN1:  gene->locus

    LocusTag:
               the LocusTag value
               ASN1:  gene->locus-tag

    Synonyms:
               bar-delimited set of unofficial symbols for the gene

    dbXrefs:
               bar-delimited set of identifiers in other databases
               for this gene.  The unit of the set is database:value.

    chromosome:
               the chromosome on which this gene is placed.
               for mitochondrial genomes, the value 'MT' is used.

    map location:
               the map location for this gene

    description:
               a descriptive name for this gene

    type of gene:
               the type assigned to the gene according to the list of options
               provided in http://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/lxr/source/src/objects/entrezgene/entrezgene.asn


    Symbol from nomenclature authority:
                when not '-', indicates that this symbol is from a
                a nomenclature authority

    Full name from nomenclature authority:
                when not '-', indicates that this full name is from a
                a nomenclature authority

    Nomenclature status:
                when not '-', indicates the status of the name from the
                nomenclature authority (O for official, I for interim)

    Other designations:
                pipe-delimited set of some alternate descriptions that
                have been assigned to a GeneID
                '-' indicates none is being reported.

    Modification date:
                the last date a gene record was updated, in YYYYMMDD format
    """
    reader = csv.reader(f, delimiter='\t')
    for l in reader:
        taxid = l[0]
        if taxid_filter and taxid_filter != taxid:
            continue
        symbol = l[10]
        yield Gene(
          id=int(l[1]),
          symbol=l[10],
          refs = l[5].split('|')
        )

def gene_map(genes):
    "Transforms a sequence of genes into a dict from ids to genes"
    return dict((gene.id, gene) for gene in genes)

gene_info_file = os.path.join(biopsy.get_data_dir(), 'NCBI', '/gene_info.gz')
def _mouse_genes():
    "Return a mapping from gene ids to mouse genes"
    print 'Parsing entrez gene file for mouse genes: %s' % gene_info_file
    return gene_map(parse_gene_info(gzip.open(gene_info_file)))

mouse_genes = lazy.PersistableLazyInitialiser(
        _mouse_genes,
        os.path.join(biopsy.get_data_dir(), 'identifiers', 'entrez', 'mouse_genes.pickle')
)


def parse_gene_protein_accessions(f, taxid_filter='10090'):
    """
            Parse a gene info file object. Yields tuples:

                    id, protein acc, protein gi

            E.g. ftp://ftp.ncbi.nih.gov/gene/DATA/gene2accesion.gz

    ===========================================================================
    gene2accession
    ---------------------------------------------------------------------------
               This file is a comprehensive report of the accessions that are
                 related to a GeneID.  It includes sequences from the international
                 sequence collaboration, Swiss-Prot, and RefSeq.

               This file can be considered as the logical equivalent of

                        ftp://ftp.ncbi.nih.gov/refseq/LocusLink/loc2ref
                   AND
                        ftp://ftp.ncbi.nih.gov/refseq/LocusLink/loc2acc

               tab-delimited
               one line per genomic/RNA/protein set of sequence accessions
               Column header line is the first line in the file.
    ---------------------------------------------------------------------------

    tax_id:
               the unique identifier provided by NCBI Taxonomy
               for the species or strain/isolate

    GeneID:
               the unique identifier for a gene
               --note:  for genomes previously available from LocusLink,
                        the identifiers are equivalent

    status:
                status of the RefSeq if a refseq, else '-'

    RNA nucleotide accession.version:
               may be null (-) for some genomes

    RNA nucleotide gi:
               the gi for an RNA nucleotide accession, '-' if not applicable

    protein accession.version:
               will be null (-) for RNA-coding genes

    protein gi:
               the gi for a protein accession, '-' if not applicable

    genomic nucleotide accession.version:
               may be null (-)

    genomic nucleotide gi:
               the gi for a genomic nucleotide accession, '-' if not applicable

    start position on the genomic accession:
                position of the gene feature on the genomic accession,
                '-' if not applicable
                position 0-based
                NOTE: this file not report the position of each exon
                for positions on RefSeq contigs and chromosomes,
                use the seq_gene.md file in the desired build directory.
                For example, for human at the time this was written:
                /am/ftp-genomes/H_sapiens/maps/mapview/BUILD.35.1
                WARNING: positions in these files are one-based, not
                0-based

    end position on the genomic accession:
                position of the gene feature on the genomic accession,
                '-' if not applicable
                position 0-based

                NOTE: this file not report the position of each exon
                for positions on RefSeq contigs and chromosomes,
                use the seq_gene.md file in the desired build directory.
                For example, for human at the time this was written:
                /am/ftp-genomes/H_sapiens/maps/mapview/BUILD.35.1
                WARNING: positions in these files are one-based, not
                0-based

    orientation:
                orientation of the gene feature on the genomic accession,
                '?' if not applicable

    assembly:
                the name of the assembly
                '-' if not applicable
"""
    reader = csv.reader(f, delimiter='\t')
    for l in reader:
        taxid = l[0]
        if taxid_filter and taxid_filter != taxid:
            continue
        acc = l[5]
        id = l[6]
        if acc != '-' and id != '-':
            yield int(l[1]), acc, int(id)

def build_protein_maps(accessions):
    protein_acc_map = cookbook.DictOfSets()
    gene_protein_map = cookbook.DictOfSets()
    for id, protein_acc, protein_id in accessions:
        protein_acc_map[protein_acc.split('.')[0]].add(protein_id)
        gene_protein_map[id].add(protein_id)
    return protein_acc_map, gene_protein_map

Proteins = cookbook.NamedTuple('Proteins', 'acc2id geneid2protein')

gene_accession_file = os.path.join(biopsy.get_data_dir(), 'NCBI', 'gene2accession.gz')

def _mouse_proteins():
    print 'Parsing entrez gene file for mouse accessions: %s' % gene_info_file
    acc2id, geneid2protein = build_protein_maps(
            parse_gene_protein_accessions(
                    gzip.open(gene_accession_file)))
    return Proteins(acc2id=acc2id, geneid2protein=geneid2protein)

mouse_proteins = lazy.PersistableLazyInitialiser(
        _mouse_proteins,
        os.path.join(biopsy.get_data_dir(), 'identifiers', 'entrez', 'mouse_proteins.pickle')
)
