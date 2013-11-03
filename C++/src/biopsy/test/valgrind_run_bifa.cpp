/**
 * Copyright John Reid 2013
 *
 * @file Code to test whether there is a memory leak when scoring reverse complements with
 * the BiFA algorithm.
 */

#include <biopsy/init.h>
#include <biopsy/analyse.h>
#include <biopsy/pssm.h>

#include <bio/fasta.h>

int
main( int argc, char * argv [] ) {
    using namespace biopsy;

    init();

    //
    // Limit number of boxes in maximal chain analysis
    //
    pssm_parameters::singleton().max_chain_num_boxes_limit = 100;

    //
    // Open FASTA
    //
    const std::string fasta_filename =
        argc > 1
        ? argv[ 1 ]
        : "/home/john/Dev/Bio/Biopsy/etc/test-data/maximal-chain-e2f.fasta";
    std::ifstream stream( fasta_filename.c_str() );
    if( ! stream.good() ) {
        throw std::logic_error( "Could not open FASTA file." );
    }

    //
    // Parse FASTA
    //
    BIO_NS::fasta_file_map_t fasta_seqs;
    const std::string centre_key = BIO_NS::parse_fasta_2( stream, fasta_seqs );

    //
    // Build sequence vector
    //
    sequence_vec_ptr seqs( new sequence_vec );
    BOOST_FOREACH( BIO_NS::fasta_file_map_t::value_type & x, fasta_seqs ) {
        x.second.resize( 100 ); // reduce length of sequence
        if( x.first == centre_key ) {
            seqs->insert( seqs->begin(), x.second );
        } else {
            seqs->push_back( x.second );
        }
    }
    BOOST_FOREACH( const sequence & seq, *seqs ) {
        std::cout << "Sequence has " << seq.size() << " base pairs.\n";
    }

    //
    // Use PSSM V$E2F_Q2
    //
    string_vec_ptr pssm_names( new string_vec );
    pssm_names->push_back( "M00803" ); // V$E2F_Q2

    //
    // score sequence
    //
    binding_hit::vec_ptr adjusted_hits;
    binding_hit::vec_ptr maximal_chain;
    binding_hits_vec_ptr unadjusted_hits;
    tie( adjusted_hits, maximal_chain, unadjusted_hits ) =
        score_pssms_on_phylo_sequences( pssm_names, seqs, .016 , .016 );
    BOOST_FOREACH( const binding_hit & hit, *adjusted_hits ) {
        std::cout << hit << "\n";
    }

    return 0;
}




