/* Copyright John Reid 2007
*/

#include "bio-pch.h"


#include "bio/defs.h"


#include "bio/biobase_db.h"
#include "bio/biobase_data_traits.h"
USING_BIO_NS

using namespace boost;
using namespace boost::archive;

using namespace std;



BIO_NS_START


namespace detail {

bool
never_parse( Database type )
{
	switch( type )
	{
	case AFFY_PROBE_DB:
	case BKL_DB:
	case RSNP_DB:
	case UNKNOWN_DB:
		return true;
	default:
		return false;
	}
}

} //namespace detail



AlignDesc::AlignDesc()
: positive_orientation( true )
{
}


bool
Site::parse_db_type( Database type )
{
	if( detail::never_parse( type ) ) { return false; }
	return false;
}



bool
Gene::parse_db_type( Database type )
{
	if( detail::never_parse( type ) ) { return false; }
	return true;
}


bool
Factor::parse_db_type( Database type )
{
	if( detail::never_parse( type ) ) { return false; }
	return true;
}


bool
Evidence::parse_db_type( Database type )
{
	if( detail::never_parse( type ) ) { return false; }
	return true;
}



int Site::get_transfac_gene_accession() const
{
	static const std::string tag = "Gene:";
	static const std::string dot = ".";

	//find the gene tag in description
	std::string::size_type pos = description.find( tag );
	if ( std::string::npos == pos )
	{
		return -1;
	}
	pos += tag.size();

	//ignore white space
	while ( ' ' == description[ pos ] )
	{
		++pos;
	}

	//the dot at the end of the text of the link to the gene table
	std::string::size_type end = description.find( dot, pos );

	//the text of the link to the gene table
	const std::string gene_link = description.substr( pos, end - pos );

	//parse the text
	const TableLink link = parse_table_link_accession_number( gene_link );
	if ( GENE_DATA != link.table_id )
	{
		return -1;
	}

	return link.entry_idx;
}





BIO_NS_END
