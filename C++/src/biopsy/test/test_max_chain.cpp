/**
 * Copyright John Reid 2013
 *
 * @file Code to test BiFa analysis.
 */

#include <biopsy/init.h>
#include <biopsy/analyse.h>
#include <biopsy/pssm.h>

#include <bio/fasta.h>

int
main( int argc, char * argv [] ) {
    using namespace biopsy;

    init();
    pssm_parameters::singleton().use_score = false;

    string_vec_ptr pssm_names( new string_vec );
    pssm_names->push_back( "M00803" ); // V$E2F_Q2

    //
    // Open FASTA
    //
    if( argc < 2 ) {
        throw std::invalid_argument( "USAGE: test_max_chain <FASTA file>" );
    }
    const std::string fasta_filename = argv[ 1 ];
    //const std::string fasta_filename = "/home/john/Dev/Bio/Biopsy/etc/test-data/maximal-chain-e2f.fasta";
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
        if( x.first == centre_key ) {
            seqs->insert( seqs->begin(), x.second );
        } else {
            seqs->push_back( x.second );
        }
    }
//    seqs->push_back( sequence( "GGCGCG" ) ); // consensus sequence = GGCGCG
//    seqs->push_back( sequence( "GGCGCG" ) );
//    seqs->push_back( sequence( "GGCGCG" ) );
//    seqs->push_back( sequence( "GGCGCG" ) );
    BOOST_FOREACH( const sequence & seq, *seqs ) {
        std::cout << "Sequence has " << seq.size() << " base pairs.\n";
    }

    //
    // reset max boxes
    //
    //pssm_parameters::singleton().max_chain_num_boxes_limit = 100000; // takes about 1 min 45 secs using about 91K boxes
    //pssm_parameters::singleton().max_chain_num_boxes_limit = 50000; // takes about 21 secs about 41K boxes
    pssm_parameters::singleton().max_chain_num_boxes_limit = 25000; // takes about 8 secs about 21K boxes

    //
    // score PSSMs and get maximal chain
    //
    binding_hit::vec_ptr hits( new binding_hit::vec() );
    const double threshold = .016;
    const double phylo_threshold = threshold;
    std::cout << "Running maximal chain algorithm.\n";
    binding_hit::vec_ptr adjusted_hits;
    binding_hit::vec_ptr maximal_chain;
    binding_hits_vec_ptr unadjusted_hits;
    tie( adjusted_hits, maximal_chain, unadjusted_hits ) =
        score_pssms_on_phylo_sequences( pssm_names, seqs, threshold, phylo_threshold );
    std::cout << "Got " << adjusted_hits->size() << " adjusted hits.\n";
    std::cout << "Hits over all sequences:\n";
    BOOST_FOREACH( binding_hit::vec_ptr hits, *unadjusted_hits ) {
        std::cout << "\tGot " << hits->size() << " unadjusted hits.\n";
    }
    if( maximal_chain ) {
        std::cout << "Got " << maximal_chain->size() << " hits in maximal chain.\n";
        if ( ! maximal_chain->size() ) {
            throw std::logic_error( "Did not get any hits in maximal chain.\n" );
        }
    } else {
        throw std::logic_error( "Did not get a maximal chain.\n" );
    }

    namespace fs = boost::filesystem;
    biopsy::build_svg(
        fs::path( "max-chain.svg" ),
        "Maximal chain test",
        ( *seqs )[ 0 ],
        threshold, // min_threshold
        *adjusted_hits,
        10, // max_num_factors
        false,
        false,
        maximal_chain.get(),
        "", // notes
        1. // max_threshold
    );

    return 0;
}
