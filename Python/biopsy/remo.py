#
# Copyright John Reid 2006
#

from _biopsy import *



def _aligned_sequences_str( self ):
    """Returns the aligned sequences as a string"""
    return "\n".join( str(s) for s in self.get_sequence_ids() )
Remome.AlignedSequenceSet.__str__ = _aligned_sequences_str


class RemomeAnalysis:
    def __init__( self, remome ):
        self.remome = remome
        self.aligned_seqs = self.remome.get_aligned_sequences()
        self.analysis = Analysis()

    def analyse( self, pssms, threshold, phylo_threshold ):
        """Analyses all those remos that aren't already in the analysis"""
        for al in self.aligned_seqs:
            for r in self.remome.get_remos_for( al ):
                k = al.get_remo_id( r )
                try:
                    hits = self.analysis.get_hits_for( k )
                except:
                    seqs = al.get_sequences_for_remo( r )
                    print 'Analysing; # seqs: %d; %s' % ( len(seqs), k )
                    hits = score_pssms_on_phylo_sequences(
                            pssms,
                            seqs,
                            threshold,
                            phylo_threshold
                    )[0]
                    self.analysis.set_hits_for( k, hits )

    def recover_from_disk(
            filename
    ):
        # try and load any previous results
        try:
            print 'Recovering from analysis from %s' % ( filename )
            self.analysis = Analysis.deserialise( filename )
            print 'Loaded previous analysis from %s' % ( filename )
            print 'Have analysis for %d remos' % ( len( self.analysis.get_keys() ) )
        except:
            pass

    def analyse_and_serialise(
            self,
            filename,
            pssms,
            threshold,
            phylo_threshold
    ):
        """Analyses all the remomes and serialises results to filename"""
        try:
            self.analyse(
                    pssms,
                    threshold,
                    phylo_threshold
            )
        finally:
            print 'Saving to %s' % ( filename )
            self.analysis.serialise( filename )
            print 'Have analysis for %d remos' % ( len( self.analysis.get_keys() ) )


def check_remome_uniqueness( r ):
    """Checks to see which sequences in a remome are identical"""
    hashed_sequences = { }
    num_same_hash = 0
    num_identical = 0
    num_false_alarms = 0
    for s in r.get_aligned_sequences():
        for remo in r.get_remos_for( s ):
            seq = remo.get_sequence_for( s.centre_sequence, True )
            h = seq.__hash__()
            if hashed_sequences.has_key( h ):
                num_same_hash += 1
                other_s = hashed_sequences[ h ][ 0 ]
                other_remo = hashed_sequences[ h ][ 1 ]
                print '%s:%s:%s and\n%s:%s:%s have the same hash\n' % (
                        s.centre_sequence,
                        s.get_sequence_info( s.centre_sequence ).region,
                        remo.get_sequences( s.centre_sequence )[0].location,
                        other_s.centre_sequence,
                        other_s.get_sequence_info( other_s.centre_sequence ).region,
                        other_remo.get_sequences( other_s.centre_sequence )[0].location,
                )
                other_seq = other_remo.get_sequence_for( other_s.centre_sequence, True )
                if other_seq == seq:
                    # print 'Found a true identical match over %d bases' % ( len( seq ) )
                    num_identical += 1
                    print seq
                else:
                    # print 'False alarm'
                    num_false_alarms += 1
                if str(s.centre_sequence) == str(other_s.centre_sequence):
                    print seq
                    print other_seq
            else:
                hashed_sequences[ h ] = ( s, remo )
    print '# with same hash: %d' % ( num_same_hash )
    print '# identical: %d' % ( num_identical )
    print '# false alarms: %d' % ( num_false_alarms )

def blast_remos( r, db = 'nr' ):
    """Uses blast to find remos in a genome"""
    from Bio.Blast import NCBIWWW, NCBIXML
    import cStringIO
    b_parser = NCBIXML.BlastParser()
    E_VALUE_THRESH = 0.04
    for s in r.get_aligned_sequences():
        for remo in r.get_remos_for( s ):
            seq = remo.get_sequence_for( s.centre_sequence, False )
            print 'Blasting: %s...' % ( seq[:60] )
            result_handle = NCBIWWW.qblast( 'blastn', db, seq )
            blast_results = result_handle.read()
            blast_out = cStringIO.StringIO(blast_results)
            b_record = b_parser.parse(blast_out)
            for alignment in b_record.alignments:
                for hsp in alignment.hsps:
                    if hsp.expect < E_VALUE_THRESH:
                        print '****Alignment****'
                        print 'sequence:', alignment.title
                        print 'length:', alignment.length
                        print 'e value:', hsp.expect
                        print 'sbjct_start:', hsp.sbjct_start
                        print hsp.query[0:75] + '...'
                        print hsp.match[0:75] + '...'
                        print hsp.sbjct[0:75] + '...'
            break
        break


def create_older_gene_set( filename = 'older_genes.txt' ):
    older_genes = set()
    for line in open( filename ):
        line = line.rstrip( '\r\n' )
        ( gene, db_name ) = line.split( ' ' )
        older_genes.add( ( gene, db_name ) )
    return older_genes

def remove_older_aligned_sequences( r, older_genes ):
    print 'Removing older aligned remos'
    removed = 0
    for s in r.get_aligned_sequences():
        gene = (
                str( s.centre_sequence.gene_id ),
                str( s.centre_sequence.db_id )
        )
        # print gene
        if gene in older_genes:
            r.remove_remos_for( s )
            removed += 1
    print '# removed by older: %d' % ( removed )

def remove_older_versioned_sequences( r ):
    print 'Removing older versioned remos'
    aligned_seqs = r.get_aligned_sequences()
    to_remove = set()
    for s1 in aligned_seqs:
        for s2 in aligned_seqs:
            if s2 == s1: break
            if (
                    s2.centre_sequence.gene_id == s1.centre_sequence.gene_id
                    and s2.centre_sequence.transcript_id == s1.centre_sequence.transcript_id
                    and s2.centre_sequence.db_id == s1.centre_sequence.db_id
                    and s2.get_sequence_info( s2.centre_sequence ).region
                            == s1.get_sequence_info( s1.centre_sequence ).region
                    and s2.centre_sequence.version < s1.centre_sequence.version
            ):
                to_remove.add( s2 )
    for s in to_remove:
        r.remove_remos_for( s )
    print '# removed by version: %d' % ( len( to_remove ) )


def get_remome_centre_sequence_genes( r ):
    """Returns a set of ensembl ids"""
    return set( s.centre_sequence.gene_id for s in r.get_aligned_sequences() )
