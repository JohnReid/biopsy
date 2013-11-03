

#ifndef BIO_BIOBASE_FILTER_H
#define BIO_BIOBASE_FILTER_H

#include "bio/defs.h"
#include "bio/matrix.h"
#include "bio/site.h"

#include <boost/program_options.hpp>
#include <boost/iterator.hpp>


BIO_NS_START




/** Filters pssms in biobase. */
struct BiobasePssmFilter
{
	bool use_consensus_sequences;
	std::string species_filter;
	boost::regex name_regex;

	BiobasePssmFilter(
		bool use_consensus_sequences = true,
		const std::string & species_filter = "V",
		const std::string & name_regex_pattern = "." );

	bool operator()( Matrix::map_t::value_type matrix ) const;
	bool operator()( Site::map_t::value_type site ) const;

	static BiobasePssmFilter get_all_pssms_filter();

	struct arg_t
	{
		bool use_consensus_sequences;
		std::string species_filter;
		std::string name_regex_pattern;

		arg_t();

		void add_options( boost::program_options::options_description & options );

		BiobasePssmFilter operator()() const;
	};
};



/** An iterator that filters out matrices in Biobase we don't want to consider. */
typedef boost::filter_iterator< BiobasePssmFilter, Matrix::map_t::const_iterator > matrix_filter_it;

/** An iterator that filters out sites in Biobase we don't want to consider. */
typedef boost::filter_iterator< BiobasePssmFilter, Site::map_t::const_iterator > site_filter_it;



matrix_filter_it get_matrices_begin( const BiobasePssmFilter & filter = BiobasePssmFilter() );
matrix_filter_it get_matrices_end( const BiobasePssmFilter & filter = BiobasePssmFilter() );
site_filter_it get_sites_begin( const BiobasePssmFilter & filter = BiobasePssmFilter() );
site_filter_it get_sites_end( const BiobasePssmFilter & filter = BiobasePssmFilter() );


/** A range of matrix iterators. */
typedef boost::iterator_range< matrix_filter_it > matrix_filter_range;

/** A range of site iterators. */
typedef boost::iterator_range< site_filter_it > site_filter_range;

matrix_filter_range get_matrices( const BiobasePssmFilter & filter = BiobasePssmFilter() );
site_filter_range get_sites( const BiobasePssmFilter & filter = BiobasePssmFilter() );


BIO_NS_END



#endif //BIO_BIOBASE_FILTER_H
