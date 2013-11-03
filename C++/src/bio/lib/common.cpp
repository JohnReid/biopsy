/* Copyright John Reid 2007
*/

#include "bio-pch.h"


#include "bio/defs.h"



#include "bio/common.h"
#include "bio/environment.h"
#include "bio/factor.h"
#include "bio/biobase_db.h"
#include "bio/biobase_data_traits.h"
#include "bio/useradmin.h"


#include <string>
using namespace std;



BIO_NS_START

Database database_from_trans_data_type( TransData table )
{
	switch( table )
	{
	//Transfac
	case CELL_DATA:
	case MATRIX_DATA:
	case SITE_DATA:
	case FACTOR_DATA:
	case FRAGMENT_DATA:
	case GENE_DATA:
	case REFERENCE_DATA:
		return TRANSFAC_DB;

	//Transcompel
	case COMPEL_DATA:
	case EVIDENCE_DATA:
		return TRANSCOMPEL_DB;

	//Transpath
	case PATHWAY_DATA:
	case MOLECULE_DATA:
	case REACTION_DATA:
		return TRANSPATH_DB;

	//unknown
	case S_DATA:
	case UNKNOWN_DATA:
		return UNKNOWN_DB;
	}

	throw std::logic_error( BIO_MAKE_STRING( "Unknown enumeration value for TransData table type: " << int( table ) ) );
}

TransData parse_biobase_table(const std::string & table_name)
{
	if ("" == table_name || "Cell" == table_name) {
		return CELL_DATA;
	} else if ("C" == table_name || "Compel" == table_name) {
		return COMPEL_DATA;
	} else if ("ev" == table_name || "Evidence" == table_name) {
		return EVIDENCE_DATA;
	} else if ("T" == table_name || "Factor" == table_name) {
		return FACTOR_DATA;
	} else if ("FR" == table_name || "Fragment" == table_name) {
		return FRAGMENT_DATA;
	} else if ("G" == table_name || "Gene" == table_name) {
		return GENE_DATA;
	} else if ("M" == table_name || "Matrix" == table_name) {
		return MATRIX_DATA;
	} else if ("MO" == table_name || "Molecule" == table_name) {
		return MOLECULE_DATA;
	} else if ("CH" == table_name || "Pathway" == table_name) {
		return PATHWAY_DATA;
	} else if ("XN" == table_name || "Reaction" == table_name) {
		return REACTION_DATA;
	} else if ("RE" == table_name || "Reference" == table_name) {
		return REFERENCE_DATA;
	} else if ("s" == table_name || "S" == table_name) {
		return S_DATA;
	} else if ("R" == table_name || "Site" == table_name) {
		return SITE_DATA;
	} else {
		throw std::logic_error( BIO_MAKE_STRING( "parse_biobase_table: Unknown biobase table: " << table_name ) );
	}
}

const char * biobase_table_to_string(TransData table)
{
	switch (table) {
		case CELL_DATA: return "";
		case COMPEL_DATA: return "C";
		case EVIDENCE_DATA: return "ev";
		case FACTOR_DATA: return "T";
		case FRAGMENT_DATA: return "FR";
		case GENE_DATA: return "G";
		case MATRIX_DATA: return "M";
		case MOLECULE_DATA: return "MO";
		case PATHWAY_DATA: return "CH";
		case REACTION_DATA: return "XN";
		case REFERENCE_DATA: return "RE";
		case S_DATA: return "s";
		case SITE_DATA: return "R";
		default:
			throw std::logic_error( BIO_MAKE_STRING( "biobase_table_to_string: Unknown biobase table" << (int) table ) );
	}
}

std::string
TableLink::get_url() const
{
	stringstream stream;
	stream << BioEnvironment::singleton().biobase_url_prefix;

	switch(table_id)
	{
	case COMPEL_DATA:
	case EVIDENCE_DATA:
		stream << "transcompel/current/bin/getTRANSCompel.cgi?" << *this;
		break;

	case CELL_DATA:
	case FACTOR_DATA:
	case FRAGMENT_DATA:
	case GENE_DATA:
	case MATRIX_DATA:
	case REFERENCE_DATA:
	case SITE_DATA:
		stream << "get.cgi?" << *this;
		break;

	case REACTION_DATA:
	case MOLECULE_DATA:
	case PATHWAY_DATA:
		stream << "transpath/current/bin/get.cgi?" << *this;
		break;

	default:
		throw std::logic_error( BIO_MAKE_STRING( "Cannot generate url for this type of data: " << table_id ) );
	}

	return stream.str();
}

TableLink::TableLink(const std::string & accession_number)
{
	//find the first digit
	size_t i;
	for (i = 0; i < accession_number.size(); ++i)
	{
		if ('0' <= accession_number[i] && accession_number[i] <= '9') {
			break;
		}
	}

	table_id = parse_biobase_table(accession_number.substr(0, i));
	entry_idx = atoi(accession_number.substr(i).c_str());
}

TableLink::TableLink(
	TransData table_id,
	int entry_idx)
	: table_id(table_id)
	, entry_idx(entry_idx)
{ }

std::string 
TableLink::get_text() const
{
	std::stringstream stream;
	stream << table_id << entry_idx;
	return stream.str();
}

bool 
TableLink::operator==(const TableLink & rhs) const
{
	return table_id == rhs.table_id && entry_idx == rhs.entry_idx;
}

bool 
TableLink::operator<(const TableLink & rhs) const
{
	if (table_id < rhs.table_id) return true;
	if (rhs.table_id < table_id) return false;
	if (table_id == rhs.table_id && entry_idx < rhs.entry_idx) return true;
	return false;
}

TableLink parse_table_link_accession_number(const std::string & accession_number)
{
	return TableLink( accession_number );
}

std::string
FactorLink::get_text() const {
	std::stringstream stream;
	stream << link.get_text() << " " << name << " ";
	BOOST_FOREACH( const std::string & aspecies, species )
		{stream << aspecies << " ";}
	return stream.str();
}

bool
FactorLink::operator==( const FactorLink & rhs ) const 
{
	if( link != rhs.link ) return false;
	return true;
}

bool
FactorLink::operator<( const FactorLink & rhs ) const 
{
	if( link < rhs.link ) return true;
	return false;
}

const std::string Identifier::get_text() const 
{
	return species_group + "$" + factor + "_" + discriminating_extension; 
}


std::ostream &
operator<<( std::ostream & os, const TableLink & tl )
{
	switch (tl.table_id)
	{
		case CELL_DATA: os << "Cell" << setw(4); break;
		case COMPEL_DATA: os << "C" << setw(5); break;
		case EVIDENCE_DATA: os << "ev" << setw(5); break;
		case FACTOR_DATA: os << "T" << setw(5); break;
		case FRAGMENT_DATA: os << "FR" << setw(7); break;
		case GENE_DATA: os << "G" << setw(6); break;
		case MATRIX_DATA: os << "M" << setw(5); break;
		case MOLECULE_DATA: os << "MO" << setw(9); break;
		case PATHWAY_DATA: os << "CH" << setw(9); break;
		case REACTION_DATA: os << "XN" << setw(9); break;
		case REFERENCE_DATA: os << "RE" << setw(7); break;
		case S_DATA: os << "<unknown: s>" << setw(9); break;
		case SITE_DATA: os << "R" << setw(5); break;
		case UNKNOWN_DATA: os << "<unknown>" << setw(9); break;
		default:
			throw std::logic_error( BIO_MAKE_STRING( "Unknown biobase table: " << int( tl.table_id ) ) );
	}

	const char prev_fill = os.fill('0');
	os << tl.entry_idx;
	os.fill(prev_fill);

	return os;
}

bool
BiobaseTablePssmEntry::is_vertebrate() const
{
	for (FactorLinkList::const_iterator f = get_factors().begin();
		get_factors().end() != f;
		++f)
	{
		Factor * factor = BiobaseDb::singleton().get_entry< FACTOR_DATA >(f->get()->link);
		if (0 != factor && factor->is_vertebrate())
		{
			return true;
		}
	}

	return false;
}

std::string TableLink::get_name() const
{
	return BiobaseDb::singleton().get_entry(*this)->get_name();
}


bool contains_factor( const FactorLinkList & factor_links, TableLink factor )
{
	for( FactorLinkList::const_iterator f = factor_links.begin();
		factor_links.end() != f;
		++f )
	{
		if( f->get()->link == factor )
		{
			return true;
		}
	}

	return false;
}

BIO_NS_END
