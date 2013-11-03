/* Copyright John Reid 2007
*/

#include "bio-pch.h"


#include "bio/defs.h"
#include "bio/factor.h"
#include "bio/biobase_db.h"
#include "bio/biobase_data_traits.h"

#include <boost/algorithm/string/predicate.hpp>

BIO_NS_START


Factor::type parse_factor_type( const std::string & s )
{
	using namespace boost::algorithm;
    if( iequals( s, "family" ) ) return Factor::family_type;
    if( iequals( s, "isogroup" ) ) return Factor::isogroup_type;
    if( iequals( s, "basic" ) ) return Factor::basic_type;
    if( iequals( s, "complex" ) ) return Factor::complex_type;
    if( iequals( s, "miRNA basic" ) ) return Factor::miRNA_type;
    if( iequals( s, "family_mod" ) ) return Factor::family_mod_type;
    if( iequals( s, "isogroup_mod" ) ) return Factor::isogroup_mod_type;
    if( iequals( s, "basic_mod" ) ) return Factor::basic_mod_type;
    if( iequals( s, "complex_mod" ) ) return Factor::complex_mod_type;
    if( iequals( s, "pre-miRNA basic" ) ) return Factor::pre_miRNA_type;
	return Factor::unknown_type;
}


const TableLink & BiobaseFactorKeyTransformer::operator()( const TableLink & link ) const
{
	return link;
}


const TableLink & BiobaseFactorKeyTransformer::operator()( const Factor * ptr ) const
{
	return ptr->accession_number;
}


const TableLink & BiobaseFactorKeyTransformer::operator()( const Factor::ptr_t ptr ) const
{
	return ptr->accession_number;
}

/**
const TableLink & BiobaseFactorKeyTransformer::operator()( const Factor::map_t::value_type factor ) const
{
	return factor.first;
}
*/

const TableLink & BiobaseFactorKeyTransformer::operator()( const FactorLink & factor_link ) const
{
	return factor_link.link;
}


const TableLink & BiobaseFactorKeyTransformer::operator()( const FactorLinkPtr & factor_link ) const
{
	return factor_link->link;
}

bool CaseInsensitiveCmp::operator()(const TableLink & s1, const TableLink & s2) const
{
	const Factor * f1 = BiobaseDb::singleton().get_entry< FACTOR_DATA >(s1);
	const Factor * f2 = BiobaseDb::singleton().get_entry< FACTOR_DATA >(s2);

	if (0 == f1 && 0 == f2)
	{
		return false;
	}
	else if (0 == f1)
	{
		return true;
	}
	else if (0 == f2)
	{
		return false;
	}
	else
	{
		return STRICMP(f1->get_name().c_str(), f2->get_name().c_str()) < 0;
	}
}





BIO_NS_END
