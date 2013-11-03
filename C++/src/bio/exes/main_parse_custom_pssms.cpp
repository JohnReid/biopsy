/* Copyright John Reid 2006, 2007
*/

#include "bio-pch.h"


#include "bio/application.h"
#include "bio/environment.h"
#include "bio/biobase_filter.h"

#include "biopsy/custom_pssm.h"
#include "biopsy/transfac.h"
using namespace biopsy;


USING_BIO_NS;
using namespace biopsy;
using namespace boost;
using namespace std;


/**
Parses the custom pssm ids provided or all in the custom pssm directory.
*/
struct ParseCustomPssmsApp : Application
{
	vector< string > ids;
	bool parse_transfac;

	/** Get all the ids for which there are custom pssm files in the custom pssm dir. */
	void populate_ids_from_custom_pssm_dir()
	{
		pssm_set_ptr pssms = all_custom_pssms();
		std::copy( pssms->begin(), pssms->end(), back_inserter( ids ) );
	}

	void parse_all_transfac() const
	{
		string_vec_ptr pssm_names = get_transfac_pssm_accessions( BiobasePssmFilter::get_all_pssms_filter() );
		BOOST_FOREACH( const string & name, *pssm_names ) get_pssm( name );
	}

	ParseCustomPssmsApp()
		: parse_transfac( false )
	{
		namespace po = boost::program_options;

		get_options().add_options()
			( "custom-pssms", po::value( &ids ), "custom PSSM ids" )
			( "parse-transfac,t", po::value( &parse_transfac ), "preprocess TRANSFAC PSSMs" )
			;

		get_positional_options().add( "custom-pssms", -1 );
   	}

	int task()
	{
		if( ! ids.size() && ! parse_transfac ) populate_ids_from_custom_pssm_dir();
		if( parse_transfac ) parse_all_transfac();

		BOOST_FOREACH( const std::string & id, ids )
		{
			std::cout << "Parsing: " << id << "\n";
			custom_pssm::ptr custom_pssm = get_custom_pssm( id );
			pssm_info pssm = get_pssm( id );
		}

		save_pssm_cache_state();

		return 0;
	}
};




int
main(
	int argc,
	char * argv[] )
{
	return ParseCustomPssmsApp().main( argc, argv );
}



