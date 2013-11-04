#
# Copyright John Reid 2010
#

"""
Code to parse ORegAnno files.
"""


import xml.sax, os, logging, tempfile, cPickle, biopsy.dbs.das as das


oreganno_data_dir = "/home/john/Data/ORegAnno"
current_data_file = "cron.saved.1-Apr-2010.xml"
oreganno_genomes = {
    'hg18'     : "http://genome.cse.ucsc.edu:80/cgi-bin/das/hg18",
    # ('hg19',    "http://genome.cse.ucsc.edu:80/cgi-bin/das/hg19"), # hg19 has no ORegAnno track??
    'mm9'      : "http://genome.cse.ucsc.edu:80/cgi-bin/das/mm9",
    'dm2'      : "http://genome.cse.ucsc.edu:80/cgi-bin/das/dm2",
    'dm3'      : "http://genome.cse.ucsc.edu:80/cgi-bin/das/dm3",
}
"Dict mapping genomes to UCSC DAS URLs."


def full_path(filename):
    "@return: The full path for the filename."
    return os.path.join(oreganno_data_dir, filename)


class Record(object):
    "Represents a record in the ORegAnno database."
    pass


class Evidence(object):
    "Represents a piece of evidence for a record in the ORegAnno database."
    pass


def is_chip_record(record):
    "@return: True iff record is based on ChIP evidence."
    try: record.evidenceSet
    except AttributeError: return False
    for evidence in record.evidenceSet:
        if 'OREGET00003' == evidence.typeStableId:
            return True
    return False


class HandleContent(xml.sax.handler.ContentHandler):
    "Handle content from ORegAnno file."

    def __init__(self, records):
        self.records = records

    def startElement(self, name, attrs):
        # logging.debug("%s: %s", name, attrs)
        if 'record' == name:
            self.record = Record()
        if 'evidenceSet' == name:
            self.evidenceSet = set()
        if 'evidence' == name:
            self.evidence = Evidence()

    def endElement(self, name):
        if 'tfName' == name:
            self.record.tfName = self.content
        elif 'stableId' == name:
            self.record.stableId = self.content
        elif 'speciesName' == name:
            self.record.speciesName = self.content
        elif 'outcome' == name:
            self.record.outcome = self.content
        elif 'evidenceClassStableId' == name:
            self.evidence.classStableId = self.content
        elif 'evidenceSubtypeStableId' == name:
            self.evidence.subtypeStableId = self.content
        elif 'evidenceTypeStableId' == name:
            self.evidence.typeStableId = self.content
        elif 'evidence' == name:
            self.evidenceSet.add(self.evidence)
            del self.evidence
        elif 'evidenceSet' == name:
            self.record.evidenceSet = self.evidenceSet
            del self.evidenceSet
        elif 'record' == name:
            self.records.append(self.record)
            del self.record

    def characters(self, content):
        self.content = content


def parse_oreganno_xml(xml_filename=None):
    "Parse an oreganno XML file."
    if None == xml_filename:
        xml_filename = full_path(current_data_file)
    tmp_file = None
    records = []
    for l in open(xml_filename):
        if l.startswith('<?'):
            if None != tmp_file:
                tmp_file.seek(0)
                content_handler = HandleContent(records)
                xml.sax.parse(tmp_file, content_handler)
            logging.debug('Creating temporary file')
            tmp_file = tempfile.TemporaryFile('r+')
        else:
            tmp_file.write(l)
    return records


def get_ucsc_oreganno_track(url):
    "@return: Yield regions corresponding to the ORegAnno track in UCSC indexed by the URL."
    logging.info('Enumerating entry points to %s', url)
    entry_points = [id for id, start, stop, orientation in das.fetch_entry_points(url)]
    logging.info('Fetching ORegAnno features from %s', url)
    oreganno_segments = das.fetch_features(url, entry_points, types=('oreganno',))
    for segment in oreganno_segments:
        id_ = segment.get('id')
        for feature in segment.findall('FEATURE'):
            label = feature.get('label')
            start = int(feature.find('START').text)
            end = int(feature.find('END').text)
            orientation = feature.find('ORIENTATION').text
            yield label, id_, start, end, orientation


def get_dna_for_region(url, label, id_, start, end, padding_bases):
    "@return: A Bio.Seq.SeqRecord for the specified region."
    segment = '%s:%d,%d' % (id_, start-padding_bases, end+padding_bases)
    _, _, _, _, sequence = das.fetch_dna(url, segment)
    logging.debug('%s: DNA from %s for segment %s: %s' % (label, url, segment, sequence))
    return SeqRecord(Seq(sequence, DNAAlphabet()), id=label, description='%s %s' % (url, segment))


def get_ucsc_track_dna(url, padding_bases):
    "@return: A list of SeqRecords for the regions in the UCSC ORegAnno track."
    return [
        get_dna_for_region(url, label, id_, start, end, padding_bases)
        for label, id_, start, end, orientation
        in get_ucsc_oreganno_track(url)
    ]


def regions_filename(genome):
    "@return: The full path of the regions FASTA file."
    return full_path("oreganno-regions-%s.fasta" % genome)



if '__main__' == __name__:
    logging.basicConfig(level=logging.INFO)

    #records = parse_oreganno_xml()
    #cPickle.dump(records, open(os.path.join(oreganno_data_dir, 'oreganno.pickle'), 'w'), cPickle.HIGHEST_PROTOCOL)

    padding_bases = 100
    for genome, url in oreganno_genomes.iteritems():
        fasta_filename = regions_filename(genome)
        if os.path.exists(fasta_filename):
            logging.info('%s already exists, skipping', fasta_filename)
        else:
            regions = get_ucsc_track_dna(url, padding_bases)
            logging.info('Writing %d regions to %s', len(regions), fasta_filename)
            SeqIO.write(regions, open(fasta_filename, "w"), "fasta")
