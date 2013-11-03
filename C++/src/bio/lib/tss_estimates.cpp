/* Copyright John Reid 2007
*/

#include "bio-pch.h"


#include "bio/defs.h"

#include "bio/tss_estimates.h"
#include "bio/environment.h"
#include "bio/serialisable.h"

#include <boost/tokenizer.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/fstream.hpp>
namespace fs = boost::filesystem;

BIO_NS_START




TssEstimate::Exon::Exon( int start, int end )
: start( start )
, end( end )
{
}


TssEstimate & TssEstimate::parse( const std::string & str )
{
	typedef boost::tokenizer< boost::char_separator< char > > tokenizer_t;

	boost::char_separator< char > sep( ";", "", boost::keep_empty_tokens );
	tokenizer_t tokenizer( str, sep );
	


	tokenizer_t::iterator tok_iter = tokenizer.begin();
	if( tok_iter == tokenizer.end() )
	{
		throw std::logic_error( "TssEstimate::parse: Expecting id" );
	}



	if( ++tok_iter == tokenizer.end() )
	{
		throw std::logic_error( "TssEstimate::parse: Expecting gene" );
	}
	gene = *tok_iter;



	if( ++tok_iter == tokenizer.end() )
	{
		throw std::logic_error( "TssEstimate::parse: Expecting transcript" );
	}
	transcript = *tok_iter;



	if( ++tok_iter == tokenizer.end() )
	{
		throw std::logic_error( "TssEstimate::parse: Expecting database" );
	}
	database = *tok_iter;



	if( ++tok_iter == tokenizer.end() )
	{
		throw std::logic_error( "TssEstimate::parse: Expecting end" );
	}
	five_prime = ( "5prime" == *tok_iter );
	BOOST_ASSERT( "5prime" == *tok_iter || "3prime" == *tok_iter );



	if( ++tok_iter == tokenizer.end() )
	{
		throw std::logic_error( "TssEstimate::parse: Expecting status" );
	}
	status = parse_status( *tok_iter );


#if 0
	if( NO_CONCLUSION == status || NO_CLONES == status )
	{
		std::cout << str << "\n";
		return *this; //nothing more to do
	}
#endif


	if( ++tok_iter == tokenizer.end() )
	{
		throw std::logic_error( "TssEstimate::parse: Expecting coordinate system" );
	}
	coord_system = parse_coord_system( *tok_iter );



	if( ++tok_iter == tokenizer.end() )
	{
		throw std::logic_error( "TssEstimate::parse: Expecting sequence region" );
	}
	seq_region = *tok_iter;



	if( ++tok_iter == tokenizer.end() )
	{
		throw std::logic_error( "TssEstimate::parse: Expecting abs_tss_pos" );
	}
	abs_tss_pos = tok_iter->empty() ? 0 : boost::lexical_cast< int >( *tok_iter );



	if( ++tok_iter == tokenizer.end() )
	{
		throw std::logic_error( "TssEstimate::parse: Expecting strand" );
	}
	positive_strand = ( "1" == *tok_iter );
	BOOST_ASSERT( "1" == *tok_iter || "-1" == *tok_iter || "" == *tok_iter );


	if( ++tok_iter == tokenizer.end() )
	{
		throw std::logic_error( "TssEstimate::parse: Expecting rel_tss_pos" );
	}
	rel_tss_pos = tok_iter->empty() ? 0 : boost::lexical_cast< int >( *tok_iter );



	if( ++tok_iter == tokenizer.end() )
	{
		throw std::logic_error( "TssEstimate::parse: Expecting first_fifty" );
	}
	seq_1st_fifty = *tok_iter;



	if( ++tok_iter == tokenizer.end() )
	{
		throw std::logic_error( "TssEstimate::parse: Expecting clone type" );
	}
	clone_type = parse_clone_type( *tok_iter );



	if( ++tok_iter == tokenizer.end() )
	{
		throw std::logic_error( "TssEstimate::parse: Expecting clone set" );
	}
	clone_set = *tok_iter;



	if( ++tok_iter == tokenizer.end() )
	{
		throw std::logic_error( "TssEstimate::parse: Expecting clone id" );
	}
	clone_id = *tok_iter;



	if( ++tok_iter == tokenizer.end() )
	{
		throw std::logic_error( "TssEstimate::parse: Expecting gene_species" );
	}
	gene_species = *tok_iter;



	if( ++tok_iter == tokenizer.end() )
	{
		throw std::logic_error( "TssEstimate::parse: Expecting clone_species" );
	}
	clone_species = *tok_iter;



	if( ++tok_iter == tokenizer.end() )
	{
		throw std::logic_error( "TssEstimate::parse: Expecting exon_list" );
	}
	{
		boost::char_separator< char > sep( "," );
		tokenizer_t exon_tokenizer( *tok_iter, sep );

		for( tokenizer_t::iterator exon_iter = exon_tokenizer.begin();
			exon_iter != exon_tokenizer.end();
			++exon_iter )
		{
			exons.push_back( Exon().parse( *exon_iter ) );
		}
	}


	if( ++tok_iter == tokenizer.end() || "" != *tok_iter )
	{
		throw std::logic_error( "TssEstimate::parse: Expecting blank token" );
	}



	if( ++tok_iter != tokenizer.end() )
	{
		throw std::logic_error( "TssEstimate::parse: Not expecting any more tokens" );
	}


	return *this;
}


TssEstimate::CoordSystem
TssEstimate::parse_coord_system( const std::string & str )
{
	if( "chromosome" == str )
	{
		return CHROMOSOME_COORDS;
	}
	else if( "scaffold" == str )
	{
		return SCAFFOLD_COORDS;
	}
	else if( "supercontig" == str )
	{
		return SUPERCONTIG_COORDS;
	}
	else if( "" == str )
	{
		return NULL_COORDS;
	}
	else 
	{
		throw std::logic_error( BIO_MAKE_STRING( "Could not parse TssEstimate::CoordSystem from \"" << str << "\"" ) );
	}
}



TssEstimate::CloneType
TssEstimate::parse_clone_type( const std::string & str )
{
	if( "" == str )
	{
		return NULL_CLONE_TYPE;
	}
	else if( "RIKEN" == str )
	{
		return RIKEN_CLONE_TYPE;
	}
	else if( "STANDARDCLONE" == str )
	{
		return STANDARD_CLONE_TYPE;
	}
	else if( "ENSEMBLTRANSCRIPT" == str )
	{
		return ENSEMBL_TRANSCRIPT_CLONE_TYPE;
	}
	else if( "RIKEN-FRAGMENT" == str )
	{
		return RIKEN_FRAGMENT_CLONE_TYPE;
	}
	else 
	{
		throw std::logic_error( BIO_MAKE_STRING( "Could not parse TssEstimate::CloneType from \"" << str << "\"" ) );
	}
}



TssEstimate::Exon &
TssEstimate::Exon::parse( const std::string & str )
{
	typedef boost::tokenizer< boost::char_separator< char > > tokenizer_t;

	boost::char_separator< char > sep( "(|)" );
	tokenizer_t tokenizer( str, sep );

	tokenizer_t::iterator it = tokenizer.begin();

	if( tokenizer.end() == it )
	{
		throw std::logic_error( "TssEstimate::parse_exon: Expecting start" );
	}
	start = it->empty() ? 0 : boost::lexical_cast< int >( *it );


	if( tokenizer.end() == ++it )
	{
		throw std::logic_error( "TssEstimate::parse_exon: Expecting end" );
	}
	end = it->empty() ? 0 : boost::lexical_cast< int >( *it );


	if( tokenizer.end() != ++it )
	{
		throw std::logic_error( "TssEstimate::parse_exon: Not expecting any more tokens" );
	}

	return *this;
}



TssEstimate::Status
TssEstimate::parse_status( const std::string & status_str )
{
	if( "NOCLONES" == status_str )
	{
		return NO_CLONES;
	}
	else if( "NOCONCLUSION" == status_str )
	{
		return NO_CONCLUSION;
	}
	else if( "CONFIRMATION" == status_str )
	{
		return CONFIRMATION;
	}
	else if( "ALTERNATIVE" == status_str )
	{
		return ALTERNATIVE;
	}
	else if( "EXONSONLY" == status_str )
	{
		return EXONS_ONLY;
	}
	else if( "PROOF" == status_str )
	{
		return PROOF;
	}
	else 
	{
		throw std::logic_error( BIO_MAKE_STRING( "Could not parse TssEstimate::Status from \"" << status_str << "\"" ) );
	}
}




std::ostream &
operator<<( std::ostream & os, TssEstimate::Status status )
{
	switch( status )
	{
	case TssEstimate::NO_CLONES:			os << "NO_CLONES"; return os;
	case TssEstimate::NO_CONCLUSION:		os << "NO_CONCLUSION"; return os;
	case TssEstimate::CONFIRMATION:			os << "CONFIRMATION"; return os;
	case TssEstimate::ALTERNATIVE:			os << "ALTERNATIVE"; return os;
	case TssEstimate::EXONS_ONLY:			os << "EXONS_ONLY"; return os;
	case TssEstimate::PROOF:				os << "PROOF"; return os;
	default:								os << "<unknown TssEstimate::Status>"; return os;
	}
}


TssEstimate::TssEstimate()
: five_prime( true )
, status( UNKNOWN_STATUS )
, coord_system( UNKNOWN_COORDS )
, abs_tss_pos( 0 )
, positive_strand( true )
, rel_tss_pos( 0 )
, clone_type( UNKNOWN_CLONE_TYPE )
{
}




std::ostream &
operator<<( std::ostream & os, TssEstimate::CoordSystem coord_system )
{
	return os;
}

std::ostream &
operator<<( std::ostream & os, TssEstimate::CloneType clone_type )
{
	return os;
}

std::ostream &
operator<<( std::ostream & os, const TssEstimate::Exon & exon )
{
	return os;
}

std::ostream &
operator<<( std::ostream & os, const TssEstimate tss_estimate )
{
	return os;
}




void
TssEstimates::parse()
{
	const fs::path path(
		BioEnvironment::singleton().get_tss_file_new_format()
	);

	fs::ifstream stream( path );

	std::string line;
	getline( stream, line ); //ignore first three lines
	getline( stream, line ); //ignore first three lines
	getline( stream, line ); //ignore first three lines
	while ( getline(stream, line) )
	{
		estimates.insert( TssEstimate().parse( line ) );
	}
}




void
TssEstimates::init_singleton()
{
	deserialise_or_init< true >(
		*this,
		fs::path(
			BioEnvironment::singleton().get_serialised_tss_estimates_file()
		),
		boost::bind< void >(
			&TssEstimates::parse,
			_1
		) 
	);
}





BIO_NS_END
