/* Copyright John Reid 2007, 2013
*/

#include "bio-pch.h"


#undef _POSIX_C_SOURCE
#include <biopsy/defs.h>
#include <biopsy/analyse.h>
#include <biopsy/custom_pssm.h>


#include "bio/application.h"
#include "bio/biobase_filter.h"
#include "bio/svg_match.h"
#include "bio/fasta.h"
USING_BIO_NS

#include <boost/range/algorithm/copy.hpp>

using namespace biopsy;
using namespace boost;
using namespace std;
namespace po = boost::program_options;


namespace biopsy {

/** Arguments for BiFa algorithm. */
struct BiFaArgs
{
    double threshold;
    double phylo_threshold;
    std::string species_filter;
    bool use_consensus_sequences;
    std::string matrix_regex;
    std::vector< std::string > pssm_sets;
    std::vector< std::string > pssm_ids;
    bool max_chain;
    unsigned max_chain_box_limit;
    double prior_odds;
    bool use_p_value;
    bool use_score;
    bool use_cumulative_dists;
    bool avg_phylo_bayes;
    double min_related_evidence;

    void add_options( po::options_description & options );

    void set_bifa_parameters() const;

    BIO_NS::BiobasePssmFilter get_filter() const;

    pssm_set_ptr get_pssm_set() const;
};


/** Display the arguments on the stream. */
std::ostream &
operator<<(std::ostream & os, const BiFaArgs & args);


pssm_set_ptr
BiFaArgs::get_pssm_set() const
{
    pssm_set_ptr result;
    // use the PSSM ids if they have been specified
    if( ! pssm_ids.empty() ) {
        result.reset( new pssm_set );
        boost::copy( pssm_ids, std::inserter( *result, result->begin() ) );
    } else {
        result = build_pssm_set_from_names( pssm_sets, use_consensus_sequences, species_filter, matrix_regex );
    }
    return result;
}


BIO_NS::BiobasePssmFilter
BiFaArgs::get_filter() const
{
    return
        BIO_NS::BiobasePssmFilter(
            use_consensus_sequences,
            species_filter,
            matrix_regex
        );
}


void
BiFaArgs::set_bifa_parameters() const {
    pssm_parameters & params = pssm_parameters::singleton();
    params.use_p_value = use_p_value;
    params.use_score = use_score;
    params.use_cumulative_dists = use_cumulative_dists;
    params.binding_background_odds_prior = prior_odds;
    params.avg_phylo_bayes = avg_phylo_bayes;
    params.min_related_evidence_fraction = min_related_evidence;
}


void
BiFaArgs::add_options( po::options_description & options )
{
    const pssm_parameters & params = pssm_parameters::singleton();
    options.add_options()
        (
            "threshold,t",
            po::value( &threshold )->default_value( .01 ),
            "Threshold for hits"
        )
        (
            "phylo-threshold,p",
            po::value( &phylo_threshold )->default_value( .01 ),
            "Threshold for hits in phlyogenetically conserved sequences"
        )
        (
            "prior-odds",
            po::value( &prior_odds )->default_value( params.binding_background_odds_prior ),
            "Threshold for hits in phlyogenetically conserved sequences"
        )
        (
            "consensus",
            po::bool_switch( &use_consensus_sequences )->default_value( false ),
            "Use consensus_sequences"
        )
        (
            "species-filter",
            po::value( &species_filter )->default_value( "V" ),
            "Species filter"
        )
        (
            "matrix_regex,r",
            po::value( &matrix_regex )->default_value( "." ),
            "Regex to match matrix names"
        )
        (
            "pssm-set,s",
            po::value( &pssm_sets ),
            "Names of custom PSSM sets"
        )
        (
            "pssm-id,i",
            po::value( &pssm_ids ),
            "Ids of PSSMs"
        )
        (
            "max-chain,m",
            po::bool_switch( &max_chain )->default_value( true ),
            "Calculate maximal chain"
        )
        (
            "max-chain-box-limit",
            po::value( &max_chain_box_limit )->default_value( params.max_chain_num_boxes_limit ),
            "Limit on number of boxes used to calculate maximal chain"
        )
        (
            "use-p-value",
            po::bool_switch( &use_p_value )->default_value( params.use_p_value ),
            "Use p-values"
        )
        (
            "use-score",
            po::bool_switch( &use_score )->default_value( params.use_score ),
            "Use score instead of BiFA"
        )
        (
            "cumulative-dists",
            po::bool_switch( &use_cumulative_dists )->default_value( params.use_cumulative_dists ),
            "Use cumulative distributions"
        )
        (
            "avg-phylo-bayes",
            po::bool_switch( &avg_phylo_bayes )->default_value( params.avg_phylo_bayes ),
            "Average the Bayes factors when adjusting for phylogenetic sequences."
        )
        (
            "min-related-evidence",
            po::value( &min_related_evidence )->default_value( params.min_related_evidence_fraction ),
            "Fraction of evidence for binding that is used as minimum in related sequences."
        )
        ;
}


std::ostream &
operator<<(std::ostream & os, const BiFaArgs & args)
{
    os
        << "Threshold: " << args.threshold << "\n"
        << "Phylo threshold: " << args.phylo_threshold << "\n"
        << "Consensus: " << (args.use_consensus_sequences ? "yes" : "no") << "\n"
        << "Species filter: " << args.species_filter << "\n"
        << "Matrix regex: " << args.matrix_regex << "\n"
        ;

    os << "Pssm sets: ";
    std::copy( args.pssm_sets.begin(), args.pssm_sets.end(), ostream_iterator< string >( os, "," ) );
    os << "\n";

    os << "Pssm IDs: ";
    std::copy( args.pssm_ids.begin(), args.pssm_ids.end(), ostream_iterator< string >( os, "," ) );
    os << "\n";

    os
        << "Prior odds: " << args.prior_odds << "\n"
        << "Maximal chain: " << (args.max_chain ? "yes" : "no") << "\n"
        << "Maximal chain box limit: " << args.max_chain_box_limit << "\n"
        << "p-value: " << (args.use_p_value ? "yes" : "no") << "\n"
        << "Use score (instead of BiFA algorithm): " << (args.use_score ? "yes" : "no") << "\n"
        << "Use cumulative distributions: " << (args.use_cumulative_dists ? "yes" : "no") << "\n"
        << "Average Bayes factors for phylogenetic adjustment: " << (args.avg_phylo_bayes ? "yes" : "no") << "\n"
        << "Min. related evidence fraction: " << args.min_related_evidence << "\n"
        ;

    return os;
}



} //namespace biopsy



struct BiFaApp : Application
{
    BiFaArgs bifa_args;
    BuildSvgArgs build_svg_args;
    std::string fasta_filename;

    BiFaApp() {
        po::options_description pssm( "PSSM scanning options" );
        bifa_args.add_options( pssm );
        get_options().add( pssm );

        po::options_description svg( "Build SVG options" );
        build_svg_args.add_options( get_options() );
        get_options().add( svg );

        get_options().add_options()
            ("fasta_file,f", po::value( &fasta_filename ), "fasta file to read sequences from")
            ;

        get_positional_options().add( "fasta_file", -1 );
    }


    virtual ~BiFaApp() { }


    int task()
    {
        cout << "Fasta file: " << fasta_filename << "\n";
        cout << bifa_args;
        cout << build_svg_args;

        //
        // Set the BiFA parameters from the arguments
        //
        bifa_args.set_bifa_parameters();

        //
        // Open the fasta file
        //
        ifstream _input_stream( fasta_filename.c_str() );
        if( ! _input_stream ) throw logic_error( BIOPSY_MAKE_STRING( "Bad input stream: \"" << fasta_filename << "\"" ) );

        //
        // Read the fasta file
        //
        fasta_file_map_t fasta_file_map;
        string centre_sequence_key = parse_fasta_2( _input_stream, fasta_file_map );

        //
        // Get the sequences
        //
        const seq_t * match_seq = &( fasta_file_map[ centre_sequence_key ] );
        sequence_vec_ptr sequences( new sequence_vec );
        vector< string > seq_ids;
        sequences->push_back( *match_seq );
        seq_ids.push_back( centre_sequence_key );
        BOOST_FOREACH( const fasta_file_map_t::value_type & v, fasta_file_map )
        {
            if (match_seq == &(v.second)) continue; //don't redo the centre sequence
            sequences->push_back( v.second );
            seq_ids.push_back( v.first );
        }
        cout << "Sequences: \n\t";
        std::copy( seq_ids.begin(), seq_ids.end(), ostream_iterator< string >( cout, "\n\t" ) );
        cout << "\n";

        //
        // Get all the pssm names from the arguments
        //
        pssm_set_ptr _pssm_set = bifa_args.get_pssm_set();
        string_vec_ptr pssm_names( new string_vec );
        std::copy( _pssm_set->begin(), _pssm_set->end(), back_inserter( *pssm_names ) );
        cout << "Using " << pssm_names->size() << " PSSMs\n";

        //
        // Score the pssms on the sequences - i.e. run BiFa
        //
        cout << "Running BiFa...\n";
        timer _timer; //time how long it takes...
        phylo_sequences_result result = score_pssms_on_phylo_sequences(
            pssm_names,
            sequences,
            bifa_args.threshold,
            bifa_args.phylo_threshold,
            bifa_args.max_chain
        );
        cout << "Took " << _timer.elapsed() << " seconds\n";
        cout << "Got " << result.get< 0 >()->size() << " hits\n";
        if( bifa_args.max_chain && sequences->size() > 1 && ! result.get< 1 >() ) //did we get a maximal chain?
            cout << "Did not get maximal chain of hits. Too many hits or perhaps too many sequences?" << endl;

        //
        // Add additional details to notes
        //
        ostringstream notes_stream;
        notes_stream << build_svg_args.notes << "\n";
        for( unsigned i = 0; seq_ids.size() != i; ++i )
            notes_stream << (*sequences)[i].size() << "bp: " << seq_ids[i] << "\n";
        notes_stream << "Fasta file: " << fasta_filename << "\n" << bifa_args;

        //
        // Build the SVG file
        //
        cout << "Building SVG in: " << build_svg_args.file << "\n";
        biopsy::build_svg(
            build_svg_args.file,
            build_svg_args.title,
            (*sequences)[0],
            bifa_args.threshold,
            *(result.get< 0 >()),
            build_svg_args.max_num_factors,
            build_svg_args.show_labels,
            build_svg_args.open_file,
            result.get< 1 >().get(),
            notes_stream.str() );

        return 0;
    }
};

int
main( int argc, char * argv[] )
{
    return BiFaApp().main( argc, argv );
}

