#ifndef BIO_BIOBASE_PARSE_H_
#define BIO_BIOBASE_PARSE_H_

#include "bio/defs.h"
#include "bio/biobase_data_traits.h"
#include "bio/biobase_db.h"
#include "bio/lexer.h"
#include "bio/serialisable.h"

#include <antlr/TokenStreamRecognitionException.hpp>

#include "MatrixParser.hpp"
#include "SiteParser.hpp"
#include "FactorParser.hpp"
#include "FragmentParser.hpp"
#include "GeneParser.hpp"
#include "PathwayParser.hpp"
#include "MoleculeParser.hpp"
#include "CompelParser.hpp"
#include "EvidenceParser.hpp"


#include <fstream>

BIO_NS_START


template< TransData type >
void
parse(
	typename DataTraits<type>::entry_t::map_t & map,
	const std::string & biobase_file = DataTraits<type>::get_biobase_file())
{
	using namespace std;
	ifstream is(biobase_file.c_str());
	if (! is) {

		throw std::logic_error( BIO_MAKE_STRING( "Could not open file: \"" << biobase_file << "\"" ) );
	
	} else {
	
		Lexer lexer(is);

		typename DataTraits<type>::parser_t parser(lexer.get_selector());
		parser.sps.selector = &lexer.get_selector();
		map.clear();
		parser.map = &map;

		boost::timer timer;
		try
		{
			parser.table();

			std::cout << "Parsed \"" << biobase_file << "\" - " << timer.elapsed() << "s\n";
		}
		catch (ANTLR_USE_NAMESPACE(antlr)TokenStreamRecognitionException & ex)
		{
			throw std::logic_error( BIO_MAKE_STRING( "Antlr parsing error: " << biobase_file << ": " << ex.toString() ) );
		}
	}
}


template< TransData type >
typename DataTraits< type >::entry_t::map_t &
get_deserialise_or_parse(
	const typename DataTraits< type >::entry_t::map_t & map )
{
	//get rid of const-ness if deserialising or parsing
	typename DataTraits< type >::entry_t::map_t & non_const_map = const_cast< typename DataTraits< type >::entry_t::map_t & >( map );

	namespace fs = boost::filesystem;

	if ( map.empty() )
	{
		const fs::path serialised_file( DataTraits< type >::get_serialised_binary_file() );

		deserialise_or_init< true >(
			non_const_map,
			serialised_file,
			boost::bind(
				parse< type >,
				boost::ref( non_const_map ),
				DataTraits<type>::get_biobase_file()
			)
		);
	}

	return non_const_map;
}






BIO_NS_END

#endif //BIO_BIOBASE_PARSE_H_
