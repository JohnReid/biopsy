#
# Copyright John Reid 2010
#

"""
Compile a data set to test BiFa algorithms from the ORegAnno database.
"""


import biopsy.data.oreganno as O, cPickle, biopsy as B, biopsy.transfac as T, random as R, logging
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from cookbook import DictOf
from biopsy.dbs import das
from biopsy.data import ucsc



def map_names_to_factors():
    "@return: A map from factor names to TRANSFAC factors."
    result = DictOf(set)
    for f in T.Factor.all():
        result[f.name.upper()].add(f)
        for s in f.synonyms:
            result[s.upper()].add(f)
    return result




class RandomSubsequenceRetriever(object):
    """Retrieve random sequences from a collection of sequences."""

    def __init__(self, sequences):
        self.sequences = sequences

    def _choose_sequence(self):
        return self.sequences[R.randint(0, len(self.sequences)-1)]

    def __call__(self, length):
        "Retrieve a random sub-sequence."
        while True:
            seq = self._choose_sequence()
            max_pos = len(seq)-1-length
            if max_pos > 0:
                pos = R.randint(0, max_pos)
                sub_seq = seq[pos:pos+length]
                return sub_seq





class DASServerRandomSequenceRetriever(object):
    """Retrieve random sequences from UCSC DAS server."""

    def __init__(self, prefix):
        "Prefix specifies a DAS prefix for genome of interest."
        self.prefix = prefix
        self.entry_points = [ep for ep in das.fetch_entry_points(prefix)]

    def _choose_entry_point(self):
        return self.entry_points[R.randint(0, len(self.entry_points)-1)]

    def __call__(self, length):
        "Retrieve a random sequence."
        while True:
            id, start, stop, orientation = self._choose_entry_point()
            max_pos = stop-1-length
            if max_pos > start:
                pos = R.randint(start, max_pos)
                segment = das.encode_segment(id, pos, pos+length)
                id, start, stop, version, seq = das.fetch_dna(self.prefix, segment)
                if seq.count('n') < length / 4: # make sure does not have too many 'N's
                    return seq


def generate_test_cases(ignore_chip_data=True):
    "@return: Yield test cases."

    #
    # Get map from ORegAnno stable ids to ORegAnno data such as TF names.
    #
    logging.info('Loading ORegAnno data.')
    oreganno_records = cPickle.load(open(O.full_path('oreganno.pickle')))
    oreganno_data = dict((record.stableId, record) for record in oreganno_records)


    #
    # Get a map from TF names to TRANSFAC factors.
    #
    logging.info('Mapping TF names to TRANSFAC factors.')
    names_to_factors = map_names_to_factors()


    #
    # Go through ORegAnno regions looking for records and relevant TRANSFAC PSSMs.
    # Build test cases from these data.
    #
    logging.info('Creating test cases')
    test_cases = []
    for genome in O.oreganno_genomes.keys():

        #
        # Get some sequences we can use as negative examples
        #
        logging.info('Loading control sequences for genome %s', genome)
        negative_sequences = list(ucsc.upstream_5000(genome))
        random_sequence_retriever = RandomSubsequenceRetriever(negative_sequences)

        #
        # Parse ORegAnno regions from FASTA file
        #
        oreganno_regions_filename = O.regions_filename(genome)
        for seq_record in SeqIO.parse(open(oreganno_regions_filename), "fasta", generic_dna):

            #
            # Decide if we want to use this region...
            #
            stable_id = seq_record.description.split()[0]
            if stable_id in oreganno_data:
                record = oreganno_data[stable_id]
                if not ignore_chip_data or not O.is_chip_record(record):
                    tfName = record.tfName.upper()
                    if tfName in names_to_factors:
                        pssms = set()
                        for f in names_to_factors[tfName]:
                            for m in f.matrices:
                                pssms.add(str(m))
                        if pssms:
                            if len(seq_record) < 5000: # cannot find negative example to match sequences longer than this
                                logging.debug(stable_id)
                                test_cases.append((record, seq_record, pssms, random_sequence_retriever(len(seq_record))))


    #
    # Yield positive and negative test cases
    #
    for record, positive_seq, pssms, negative_seq in test_cases:
        yield (record, positive_seq, pssms), True
        yield (record, negative_seq, pssms), False





if '__main__' == __name__:
    logging.basicConfig(level=logging.INFO)

    test_cases = list(generate_test_cases())
    logging.info('Created %d test cases', len(test_cases))
