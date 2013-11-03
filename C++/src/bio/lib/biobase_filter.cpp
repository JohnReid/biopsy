/* Copyright John Reid 2007
*/

#include "bio-pch.h"


#include "bio/defs.h"

#include "bio/biobase_filter.h"
#include "bio/biobase_db.h"
#include "bio/biobase_data_traits.h"
#include "bio/biobase_match.h"


namespace po = boost::program_options;


BIO_NS_START

bool
want_to_match(Matrix * matrix) {
	return
		is_matchable(matrix)
		&& matrix->id.species_group == "V";
}

bool
want_to_match(Site * site) {
	return
		is_matchable(site)
		&& site->id.factor == "CONS";
}

matrix_filter_it
get_matrices_begin( const BiobasePssmFilter & filter )
{
	return
		matrix_filter_it(
			filter,
			BiobaseDb::singleton().get_matrices().begin(),
			BiobaseDb::singleton().get_matrices().end());
}

matrix_filter_it
get_matrices_end( const BiobasePssmFilter & filter )
{
	return
		matrix_filter_it(
			filter,
			BiobaseDb::singleton().get_matrices().end(),
			BiobaseDb::singleton().get_matrices().end());
}

site_filter_it
get_sites_begin( const BiobasePssmFilter & filter )
{
	return
		site_filter_it(
			filter,
			BiobaseDb::singleton().get_sites().begin(),
			BiobaseDb::singleton().get_sites().end());
}

site_filter_it
get_sites_end( const BiobasePssmFilter & filter )
{
	return
		site_filter_it(
			filter,
			BiobaseDb::singleton().get_sites().end(),
			BiobaseDb::singleton().get_sites().end());
}



matrix_filter_range
get_matrices( const BiobasePssmFilter & filter )
{
	return
		boost::make_iterator_range(
			get_matrices_begin( filter ),
			get_matrices_end( filter ) );
}


site_filter_range
get_sites( const BiobasePssmFilter & filter )
{
	return
		boost::make_iterator_range(
			get_sites_begin( filter ),
			get_sites_end( filter ) );
}



template< typename MatrixT >
bool
match_matrix_name_regex( const boost::regex & re, MatrixT matrix )
{
	boost::smatch what;

	//get the list of factors
	const FactorLinkList & factors = matrix.second->get_factors();
	for (FactorLinkList::const_iterator f = factors.begin();
		factors.end() != f;
		++f)
	{
		//get this factor
		Factor * factor = BiobaseDb::singleton().get_entry< FACTOR_DATA >( (*f)->link );

		if (0 != factor)
		{
			//check the factor's name
			if (boost::regex_search(bio_to_upper(factor->get_name()), what, re))
			{
				return true;
			}

#if 0
			//check all the synonyms
			for (Factor::synonym_set_t::const_iterator s = factor->synonyms.begin();
				factor->synonyms.end() != s;
				++s)
			{
				if (boost::regex_search(bio_to_upper(*s), what, re))
				{
					return true;
				}
			}
#endif //0
		}
	}

	//check the matrix name
	bool result = boost::regex_search(matrix.second->get_name(), what, re);

	return result;
}



BiobasePssmFilter::BiobasePssmFilter(
	bool use_consensus_sequences,
	const std::string & species_filter,
	const std::string & name_regex_pattern )
	: use_consensus_sequences( use_consensus_sequences )
	, species_filter( species_filter )
	, name_regex( name_regex_pattern )
{
}

bool
BiobasePssmFilter::operator()( Matrix::map_t::value_type matrix ) const
{
	return
		species_filter.find( matrix.second->id.species_group[0] ) != std::string::npos
		&& match_matrix_name_regex( name_regex, matrix );
}



bool
BiobasePssmFilter::operator()( Site::map_t::value_type site ) const
{
	return
		use_consensus_sequences
		&& "CONS" == site.second->id.factor
		&& site.second->is_vertebrate()
		&& match_matrix_name_regex( name_regex, site );
}

BiobasePssmFilter
BiobasePssmFilter::get_all_pssms_filter()
{
	return 
		BiobasePssmFilter(
			true,
			"BFINPV",
			"." );
}

BiobasePssmFilter::arg_t::arg_t()
: use_consensus_sequences( true )
, species_filter( "V" )
, name_regex_pattern( "." )
{
}

void
BiobasePssmFilter::arg_t::add_options( po::options_description & options )
{
	options.add_options()
		( "consensus", po::value( &use_consensus_sequences )->default_value( true ), "use consensus_sequences" )
		( "species_filter", po::value( &species_filter )->default_value( "V" ), "species filter" )
		( "pssms", po::value( &name_regex_pattern )->default_value( "." ), "pssm name regex" )
		;
}

BiobasePssmFilter
BiobasePssmFilter::arg_t::operator()() const
{
	return
		BiobasePssmFilter(
			use_consensus_sequences,
			species_filter,
			name_regex_pattern
		);
}


BIO_NS_END
