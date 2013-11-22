/**
@file

Copyright John Reid 2007, 2013
*/

#include "bio-pch.h"




#include "bio/application.h"
#include "bio/biobase_likelihoods.h"
USING_BIO_NS

using namespace boost;
namespace po = boost::program_options;


#include <iostream>
#include <string>
using namespace std;


/**
 * Shows the likelihoods of scores for a PSSM for different scoring methods.
 */
struct PssmLikelihoodsApp : Application
{
	typedef std::vector< std::string > string_vec_t;

	unsigned acc_id;
	bool consensus_sequence;
	std::string sequence;

	PssmLikelihoodsApp()
	{
		get_options().add_options()
			( "accession_id,a", po::value( &acc_id )->default_value( 0 ), "PSSM accession id" )
			( "consensus_sequence,c", po::bool_switch( &consensus_sequence )->default_value( false ), "consensus sequence (site table)" )
			( "sequence,s", po::value( &sequence )->default_value( "" ), "score pssm on sequence" )
			;

		//get_positional_options().add("motif", -1);
	}

	int task()
	{
		//get the link to the pssm
		if( 0 == acc_id )
		{
			throw std::logic_error( "No accession id specified" );
		}
		const TableLink link( consensus_sequence ? SITE_DATA : MATRIX_DATA, acc_id );

		//get the binding likelihoods
		const unsigned total_counts = LikelihoodsCache::singleton().get_total_counts( link );
		const BiobaseLikelihoods * binding_likelihoods = LikelihoodsCache::singleton().get_score_likelihoods( link, false, false );
		const BiobaseLikelihoods * binding_or_better_likelihoods = LikelihoodsCache::singleton().get_score_likelihoods( link, false, true );
		const BiobaseLikelihoods * background_likelihoods = LikelihoodsCache::singleton().get_score_likelihoods( link, true, false );
		const BiobaseLikelihoods * background_or_better_likelihoods = LikelihoodsCache::singleton().get_score_likelihoods( link, true, true );

		//check we got the likelihoods
		if (
            ! binding_likelihoods
            || ! binding_or_better_likelihoods
            || ! background_likelihoods
            || ! background_or_better_likelihoods
        ) {
		    std::cerr << "ERROR: Could not get likelihoods!\n";
		    return -1;
		}

		//do we want to score the pssm on the sequence?
		BIO_NS::float_t score = -1.0;
		if( "" != sequence )
		{
			Pssm pssm = make_pssm( link );
			if( pssm.size() > sequence.size() )
			{
				throw std::logic_error( "Sequence is not long enough" );
			}
			score = pssm.score( sequence.begin(), false );
			std::cout << "Score: " << score << "\n";
		}
		std::cout.setf( ios::left, ios::adjustfield );
		std::cout.setf( ios::fixed, ios::floatfield );
		std::cout
			<< setw( 12 ) << "Score"
			<< setw( 12 ) << "Binding"
			<< setw( 12 ) << "or better"
			<< setw( 12 ) << "Background"
			<< setw( 12 ) << "or better"
			<< "\n";
		for( unsigned i = 0; binding_likelihoods->size() != i; ++i )
		{
			std::cout
				<< setw(12) << double( i ) / double( binding_likelihoods->size() )
				<< setw(12) << ( *binding_likelihoods )[ i ]
				<< setw(12) << ( *binding_or_better_likelihoods )[ i ]
				<< setw(12) << ( *background_likelihoods )[ i ]
				<< setw(12) << ( *background_or_better_likelihoods )[ i ]
				<< ( unsigned( score * 100 ) == i ? " ************ " : "" ) //flag the score that the pssm produced
				<< "\n";
		}
		
		return 0;
	}
};

int
main(int argc, char * argv[])
{
	return PssmLikelihoodsApp().main( argc, argv );
}

