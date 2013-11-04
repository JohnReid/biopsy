#
# Copyright John Reid 2010
#

"""
Code to connect to DAS servers.
"""

import logging



def getText(nodelist):
    "@return: The text in a nodelist."
    rc = ""
    for node in nodelist:
        if node.nodeType == node.TEXT_NODE:
            rc = rc + node.data
    return rc


def fetch(url):
    "@return: An element tree created by parsing the result from the URL."
    from xml.etree.ElementTree import parse
    import urllib
    logging.debug('Fetching and parsing %s', url)
    return parse(urllib.urlopen(url))



def fetch_data_sources(prefix):
    "@return: Yield the data sources for the prefix as tuples (id_, version, label, url, desc)."
    tree = fetch('%s/das/dsn' % prefix)
    root = tree.getroot()
    for dsn in root.findall('DSN'):
        source = dsn.find('SOURCE')
        mapmaster = dsn.find('MAPMASTER')
        desc = None
        desc_node = dsn.find('DESCRIPTION')
        if desc_node:
            desc = desc_node.text
        yield (
            source.get('id'),
            source.get('version'),
            source.text,
            mapmaster.text,
            desc
        )


def encode_segment(id, start, stop):
    "@return: A string representing the segment for use in DAS."
    return '%s:%d,%d' % (id, start, stop)


def fetch_dna(prefix, segment):
    "@return: The DNA for the given segment in a tuple (id, start, stop, version)"
    tree = fetch('%s/dna?segment=%s' % (prefix, segment))
    root = tree.getroot()
    sequence = root.find('SEQUENCE')
    dna = sequence.find('DNA')
    return sequence.get('id'), int(sequence.get('start')), int(sequence.get('stop')), sequence.get('version'), dna.text.strip().replace('\n', '')



def fetch_entry_points(prefix):
    "@return: Yield the entry points as (id, start, stop, orientation) tuples."
    tree = fetch('%s/entry_points' % prefix)
    root = tree.getroot()
    entry_points = root.find('ENTRY_POINTS')
    for segment in entry_points.findall('SEGMENT'):
        yield (
            segment.get('id'),
            int(segment.get('start')),
            int(segment.get('stop')),
            segment.get('orientation'),
        )


def fetch_types(prefix):
    "@return: Yield the types."
    tree = fetch('%s/types' % prefix)
    root = tree.getroot()
    segment = root.find('GFF').find('SEGMENT')
    for t in segment.findall('TYPE'):
        yield (
            t.get('category'),
            t.get('id'),
        )



def fetch_features(prefix, segments, types=None, categories=None, feature_ids=None, group_ids=None, categorize=False):
    "@return: This query returns the annotations across one or more segments of sequence."
    url = '%s/features?%s%s' % (
        prefix,
        ';'.join('segment=%s' % segment for segment in segments),
        types and ''.join(';type=%s' % type for type in types) or ''
    )
    tree = fetch(url)
    root = tree.getroot()
    gff = root.find('GFF')
    segments = root.findall('SEGMENT')
    if not segments: # UCSC stores SEGMENT node under GFF contrary to standard
        segments = gff.findall('SEGMENT')
    return segments


if "__main__" == __name__:

    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    from Bio.Alphabet import DNAAlphabet
    from Bio import SeqIO

    logging.basicConfig(level=logging.DEBUG)


    def test_ucsc_mysql():
        """
        Code to test direct mySQL access to UCSC.

        Nothing more than a connection and printing available tables at the moment.
        Would require more work to investigate database schema.
        """
        import MySQLdb
        db = MySQLdb.connect(host="genome-mysql.cse.ucsc.edu", user="genomep", passwd="password")

        db.query("SHOW DATABASES")
        r = db.store_result()
        print r.fetch_row(maxrows=0)

        db.query("SHOW TABLES FROM mm9")
        r = db.store_result()
        print r.fetch_row(maxrows=0)


    server = "http://genome.ucsc.edu/cgi-bin"
    logging.info('DAS DSNs from %s:' % server)
    for id_, version, label, url, desc in fetch_data_sources(server):
        logging.info('%s %s "%s" %s "%s"' % (id_, version, label, url, desc))
