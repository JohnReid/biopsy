/**
@file

Copyright John Reid 2007, 2013
*/

#include "bio-pch.h"


#include "bio/defs.h"

#include "bio/matrix_test_data.h"
#include "bio/fragment_test_data.h"
#include "bio/site_test_data.h"
#include "bio/biobase_db.h"
#include "bio/biobase_data_traits.h"
#include "bio/biobase_tf.h"
#include "bio/biobase_binding_model.h"
#include "bio/model_2_factor.h"
#include "bio/background_models.h"
#include "bio/biobase_filter.h"

#include <boost/iterator/transform_iterator.hpp>
#include <boost/assign/list_of.hpp>

#include <numeric>
#include <iostream>

BIO_NS_START


Test::~Test() { }


SiteTestData::SiteTestData(
	const SiteTestCase & test_case)
	: BiFaTestData( test_case.centre_sequence )
	, site( BiobaseDb::singleton().get_entry< SITE_DATA >( test_case.site ) )
{
	std::copy(
		test_case.sequences.begin(),
		test_case.sequences.end(),
		std::back_inserter( input.conserved_sequences ) );
}

SiteTestData::~SiteTestData() { }

Test::ptr_t
SiteTestData::get_test_for_algorithm( BiFaAlgorithm::ptr_t algorithm )
{
	return
		Test::ptr_t(
			new BiFaFactorTest(
				*this,
				site->factor_links,
				algorithm
			)
		)
		;
}



/** Generate some test data from a matrix by the following: Take a random sequence of the given length
and replace a section of it with one of the sequences used to generate the matrix. */
MatrixTestData::MatrixTestData(
	const Matrix * matrix,
	size_t seq_length)
	: matrix(matrix)
{

	//make sure we have a binding site to choose
	if (0 == matrix->align_descs.size())
	{
		throw std::logic_error( "No binding site to choose from" );
	}

	//generate a random sequence
	input.centre_sequence.reserve( seq_length );
	gen_sequence_from_random_species_dna_hmm( input.centre_sequence, seq_length );

	bool found_binding_site = false;
	for( AlignDescList::const_iterator ad = matrix->align_descs.begin();
		matrix->align_descs.end() != ad;
		++ad )
	{
		if( seq_length < ad->get()->sequence.size() //make sure the site is not longer than the sequence
			||
			ad->get()->sequence.size() != (size_t) ad->get()->length
			||
			matrix->get_size() != ad->get()->sequence.size()
			||
			0 == ad->get()->sequence.size()
			||
			ad->get()->sequence.end()
				!=
				std::find_if(
					ad->get()->sequence.begin(),
					ad->get()->sequence.end(),
					! boost::lambda::bind< bool >( is_known_nucleotide(), boost::lambda::_1 ) ) )
		{
			//not suitable
			continue;
		}
		else
		{
			found_binding_site = true;
			binding_site = *ad;
			break;
		}
	}
	if( ! found_binding_site )
	{
		throw std::logic_error( "Could not find suitable binding site" );
	}

	//choose a position in the sequence to put it
	unsigned binding_site_idx = get_uniform_index(seq_length - binding_site->sequence.size());

	//replace the substring with the site
	input.centre_sequence.replace(
		binding_site_idx,
		binding_site->sequence.size(),
		binding_site->sequence);
}


MatrixTestData::~MatrixTestData() { }


BiFaTestData::BiFaTestData(
	const BiFaInput & input )
	: input( input )
{
}


BiFaTestData::~BiFaTestData() { }


const BiFaInput &
BiFaTestData::get_input() const
{
	return input;
}



const BiFaOutput &
BiFaTestData::get_output_for( BiFaAlgorithm::ptr_t algorithm )
{
	output_map_t::iterator o = output.find( algorithm );
	if( output.end() == o )
	{
		o = output.insert( output_map_t::value_type( algorithm, ( *algorithm )( input ) ) ).first;
	}

	return *( o->second );
}



const TF::hit_set_t &
BiFaTestData::get_factors_for( BiFaAlgorithm::ptr_t algorithm )
{
	factor_map_t::iterator f = factor_hits.find( algorithm );
	if( factor_hits.end() == f )
	{
		f = factor_hits.insert( factor_map_t::value_type( algorithm, TF::hit_set_ptr_t( new TF::hit_set_t ) ) ).first;
		model_hits_2_factor_hits( get_output_for( algorithm ).hits, *( f->second ) );
	}

	return *( f->second );
}



/** Generate some test data from a matrix by the following: Take a random sequence of the given length
and replace a section of it with one of the sequences used to generate the matrix. */
FragmentTestData::FragmentTestData(
	const Fragment * fragment)
	: BiFaTestData( fragment->sequence )
	, fragment( fragment )
{
}


FragmentTestData::~FragmentTestData() { }


MatrixTest::MatrixTest(
	MatrixTestData & data,
	BiFaAlgorithm::ptr_t algorithm )
	: data( data )
	, algorithm( algorithm )
{
}


MatrixTest::~MatrixTest() { }


Test::ptr_t
MatrixTestData::get_test_for_algorithm( BiFaAlgorithm::ptr_t algorithm )
{
	return
		Test::ptr_t(
			new MatrixTest(
				*this,
				algorithm
			)
		)
		;
}



BinaryTestResults
MatrixTest::get_results( double threshold )
{

	//all of the models the algorithm applies
	BindingModel::set_t universe;
	algorithm->fill_model_universe( universe );

	//those that the test is looking for
	const BindingModel::set_t true_positives = boost::assign::list_of
		( algorithm->get_model_for( boost::any( data.matrix->get_link() ) ) )
		;

	return
		get_results_for(
			data.get_output_for( algorithm ).hits,
			threshold,
			universe,
			true_positives );

}

BinaryTestResults
get_results_for(
	const BindingModel::hit_set_t & hits,
	double threshold,
	const BindingModel::set_t & universe,
	const BindingModel::set_t & true_positives)
{
	//check each and add to result
	BinaryTestResults result;
	for( BindingModel::set_t::const_iterator m = universe.begin();
		universe.end() != m;
		++m )
	{
		const double p_binding = get_p_binding( hits, *m );
		const bool in_test_output = p_binding > threshold;

		const bool should_be_in_test_output = true_positives.end() != true_positives.find( *m );

		//update true positives, etc...
		result( in_test_output, should_be_in_test_output );
	}

	return result;
}


Test::ptr_t
FragmentTestData::get_test_for_algorithm( BiFaAlgorithm::ptr_t algorithm )
{
	return
		Test::ptr_t(
			new BiFaFactorTest(
				*this,
				fragment->factor_links,
				algorithm
			)
		)
		;
}

BiFaFactorTest::BiFaFactorTest(
	BiFaTestData & data,
	const FactorLinkList & factors,
	BiFaAlgorithm::ptr_t algorithm )
	: data( data )
	, factors( factors )
	, algorithm( algorithm )
{
}

BiFaFactorTest::~BiFaFactorTest() { }

BinaryTestResults
BiFaFactorTest::get_results( double threshold )
{

	//all of the models the algorithm applies
	BindingModel::set_t universe;
	algorithm->fill_model_universe( universe );

	//those that the test is looking for
	BindingModel::set_t true_positives;
	{
		typedef std::set< TableLink > table_link_set;
		table_link_set matrix_links;
		for( FactorLinkList::const_iterator f = factors.begin();
			factors.end() != f;
			++f )
		{
			const Factor * factor = BiobaseDb::singleton().get_entry< FACTOR_DATA >( f->get()->link );
			std::copy( factor->matrices.begin(), factor->matrices.end(), std::inserter( matrix_links, matrix_links.begin() ) );
		}
		for( table_link_set::const_iterator l = matrix_links.begin();
			matrix_links.end() != l;
			++l )
		{
			true_positives.insert( algorithm->get_model_for( boost::any( *l ) ) );
		}
	}

	return
		get_results_for(
			data.get_output_for( algorithm ).hits,
			threshold,
			universe,
			true_positives );
}


void create_test_data_from_matrices( BiFaTestData::vec_t & test_data, size_t seq_length )
{
	const matrix_filter_it matrices_begin = get_matrices_begin();
	const matrix_filter_it matrices_end = get_matrices_end();

	for (matrix_filter_it i = matrices_begin;
		matrices_end != i;
		++i)
	{
		try
		{
			test_data.push_back( BiFaTestData::ptr_t( new MatrixTestData( i->second.get(), seq_length) ) );
		}
		catch ( const std::exception & )
		{
			//cerr << "Error: " << i->second->get_name() << ": " << ex << endl;
		}
	}
}

void create_test_data_from_fragments( BiFaTestData::vec_t & test_data )
{
	for( Fragment::map_t::const_iterator i = BiobaseDb::singleton().get_fragments().begin();
		BiobaseDb::singleton().get_fragments().end() != i;
		++i)
	{
		test_data.push_back( BiFaTestData::ptr_t( new FragmentTestData( i->second.get() ) ) );
	}
}


void 
create_test_data_from_sites( 
	BiFaTestData::vec_t & test_data )
{
	for( SiteTestCases::const_iterator s = SiteTestCases::singleton().begin();
		SiteTestCases::singleton().end() != s;
		++s )
	{
		test_data.push_back(
			BiFaTestData::ptr_t(
				new SiteTestData(
					**s
				)
			)
		)
		;
	}
}



struct threshold_calculator
{
	double threshold;
	double max_difference;

	threshold_calculator(  )
		: threshold( 0.0 )
		, max_difference( 0.0 )
	{
	}

	void operator()(
		double this_threshold,
		double this_value,
		double next_threshold,
		double next_value)
	{
		double difference = fabs( next_value - this_value );
		if( difference > max_difference )
		{
			max_difference = difference;
			threshold = ( this_threshold + next_threshold ) / 2.0;
		}
	}
};


ROCPoint
get_roc_point(
	Test::vec_t & test_data,
	double threshold )
{
	//get the roc point for this threshold
	BinaryTestResults results =
		std::accumulate(
			boost::make_transform_iterator(
				test_data.begin(),
				boost::bind(
					&Test::get_results,
					_1,
					threshold
				)
			),
			boost::make_transform_iterator(
				test_data.end(),
				boost::bind(
					&Test::get_results,
					_1,
					threshold
				)
			),
			BinaryTestResults()
		);

	return results.get_roc_point();
}




BiFaTestData2Test::BiFaTestData2Test( BiFaAlgorithm::ptr_t algorithm )
	: algorithm( algorithm )
{
}

Test::ptr_t
BiFaTestData2Test::operator()( BiFaTestData::ptr_t data ) const
{
	return data->get_test_for_algorithm( algorithm );
}


void
run_tests(
	Test::vec_t & test_data,
	test_result_map_t & result_map,
	double min_threshold,
	double max_threshold,
	unsigned num_thresholds )
{
	result_map.clear();
    result_map.insert( test_result_map_t::value_type( min_threshold, get_roc_point( test_data, min_threshold ) ) );
    result_map.insert( test_result_map_t::value_type( max_threshold, get_roc_point( test_data, max_threshold ) ) );

	bool using_spec = false;
	while( result_map.size() < num_thresholds )
	{
		//calculate new threshold
		threshold_calculator calc;
		for( test_result_map_t::const_iterator r1 = result_map.begin();
			result_map.end() != r1;
			++r1 )
		{
			test_result_map_t::const_iterator r2 = boost::next( r1 );
			if( result_map.end() == r2 )
			{
				break;
			}

			calc(
				r1->first,
				using_spec
					? r1->second.specificity
					: r1->second.sensitivity,
				r2->first,
				using_spec
					? r2->second.specificity
					: r2->second.sensitivity );

		}
		const double threshold = calc.threshold;
		using_spec = ! using_spec;

		ROCPoint roc_point = get_roc_point( test_data, threshold );
		result_map.insert( test_result_map_t::value_type( threshold, roc_point ) );
	}
}

std::ostream &
operator<<( std::ostream & os, const test_result_map_t::value_type & value )
{
	return
		os
			<< value.first 
			<< "," << value.second.sensitivity 
			<< "," << value.second.specificity
			;
}

std::ostream &
operator<<( std::ostream & os, const test_result_map_t & result_map )
{
	os << "Threshold,Sensitivity,Specificity\n";
	std::copy(
		result_map.begin(),
		result_map.end(),
		std::ostream_iterator< test_result_map_t::value_type >( os, "\n" ) );

	return os;
}

BIO_NS_END

