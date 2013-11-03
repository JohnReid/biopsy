

#include "bio/pssm_bayesian_binding_model.h"
#include "bio/biobase_binding_model.h"
#include "bio/biobase_score.h"
#include "bio/model_2_factor.h"
#include "bio/serialisable.h"
#include "bio/biobase_filter.h"
USING_BIO_NS

#include <boost/test/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/assign/list_of.hpp>
using namespace boost;
using namespace boost::assign;
using boost::unit_test::test_suite;
namespace fs = boost::filesystem;

#include <iostream>
using namespace std;


//#define VERBOSE_CHECKING


const std::vector< TableLink > pssm_links =
	list_of
		( TableLink( MATRIX_DATA, 328 ) )
		;

const SeqList sequences =
	list_of
		( "TGACTCATGCGTAGAGAT" )
		;


template< typename It >
void
check_strict_weak_ordering(
	It begin,
	It end )
{
	for( It i1 = begin; end != i1; ++i1 )
	{
		BOOST_CHECK( ! ( *i1 < *i1 ) ); //Irreflexivity

		for( It i2 = begin; end != i2; ++i2 )
		{
			if ( *i1 < *i2 )
			{
				BOOST_CHECK( ! ( *i2 < *i1 ) ); //Antisymmetry 
				
				for( It i3 = begin; end != i3; ++i3 )
				{
					if ( *i2 < *i3 )
					{
						BOOST_CHECK( *i1 < *i3 ); //Transitivity 
					}
				}
			}

		}
	}
}

template< typename Binder >
void
check_strict_weak_ordering( const typename BindingHitSet< Binder >::type & hits )
{
	check_strict_weak_ordering( hits.get< typename BindingHitSet< Binder >::position >().begin(), hits.get< 0 >().end() );
	check_strict_weak_ordering( hits.get< typename BindingHitSet< Binder >::binder >().begin(), hits.get< 1 >().end() );
	check_strict_weak_ordering( hits.get< typename BindingHitSet< Binder >::prob >().begin(), hits.get< 2 >().end() );
}


template< typename Binder >
void display_hits( const typename BindingHitSet< Binder >::type & hits )
{
	std::copy(
		hits.get< typename BindingHitSet< Binder >::binder >().begin(),
		hits.get< typename BindingHitSet< Binder >::binder >().end(),
		std::ostream_iterator< typename BindingHitSet< Binder >::type::value_type >( std::cout, "\n" ) );
	std::cout << "\n\n";
}

void
check_model_2_factor( const BindingModel::hit_set_t & hits )
{
	TF::hit_set_t factor_hits;
	model_hits_2_factor_hits( hits, factor_hits );

#ifdef VERBOSE_CHECKING
	display_hits< TF >( factor_hits );
#endif //VERBOSE_CHECKING
}

struct CheckEqual
{
	template< typename T >
	void
	operator()( const T & t1, const T & t2 ) const
	{
		BOOST_CHECK_EQUAL( t1, t2 );
	}
};


void
check_serialisation( const BindingModel::hit_set_t & hits )
{
#ifndef _DEBUG
	check_strict_weak_ordering< BindingModel >( hits );
#endif

	fs::path test_file( "hit_serialisation_test.txt" );
	serialise< BindingModel, false >( hits, test_file );

	BindingModel::hit_set_t copy_of_hits;
	deserialise< BindingModel, false >( copy_of_hits, test_file );

	BOOST_REQUIRE( copy_of_hits.size() == hits.size() );
	BOOST_CHECK( std::equal( hits.begin(), hits.end(), copy_of_hits.begin() ) );

	std::pair< BindingModel::hit_set_t::const_iterator, BindingModel::hit_set_t::const_iterator > mm =
		std::mismatch( hits.begin(), hits.end(), copy_of_hits.begin() );
}

void
check_deaf_binding_model()
{
	cout << "******* check_deaf_binding_model()\n";

	const seq_t sequence = "NNNNNNNNNNNNNNNNNNNNNNNNN";

	BindingModel * deaf_model = BiobaseBindingModel::parameter_t( TableLink( MATRIX_DATA, 1001 ) ).get_model();

	BindingModel::hit_set_t hits;

	make_sequence_scorer(
		sequence.begin(),
		sequence.end(),
		1e-5,
		std::inserter( hits, hits.begin() )
	) (
		deaf_model
	);
}

void
check_pssm_bayesian_binding_model( TableLink link )
{
	cout << "******* check_pssm_bayesian_binding_model(): " << link << "\n";

	const seq_t sequence = "TGACTCATGCGTAGAGAT";

	BindingModel::hit_set_t hits;

	make_sequence_scorer(
		sequence.begin(),
		sequence.end(),
		1e-5,
		std::inserter( hits, hits.begin() )
	) (
		BiobaseBindingModel::parameter_t( link ).get_model()
	);

#ifdef VERBOSE_CHECKING
	display_hits< BindingModel >( hits );
#endif

	check_serialisation( hits );
	check_model_2_factor( hits );
}

/** Jerome had problem with this binding site and matrix 962. */
void
check_ar_binding_model(  )
{
	cout << "******* check_ar_biobase_binding_model()\n";

	const seq_t sequence = "AGAGCATGG" ;
	const TableLink link( MATRIX_DATA, 962 );

	BindingModel::hit_set_t hits;

	BindingModel * model = BiobaseBindingModel::parameter_t( link ).get_model() ;

	make_sequence_scorer(
		sequence.begin(),
		sequence.end(),
		1e-5,
		std::inserter( hits, hits.begin() )
	) (
		model
	);

#ifdef VERBOSE_CHECKING
	display_hits< BindingModel >( hits );
#endif

}

void
check_biobase_binding_model( const seq_t & sequence )
{
	cout << "******* check_biobase_binding_model(): " << sequence << "\n";

	BindingModel::hit_set_t hits;

	score_all_biobase_pssms(
		make_sequence_scorer(
			sequence.begin(),
			sequence.end(),
			1e-5,
			std::inserter( hits, hits.begin() ) ),
		BiobasePssmFilter(),
		Link2BiobaseBindingModel() );

#ifdef VERBOSE_CHECKING
	display_hits< BindingModel >( hits );
#endif

	check_serialisation( hits );
	check_model_2_factor( hits );
}

void
check_binding_hit_less_than()
{
	cout << "******* check_binding_hit_less_than()\n";

	BindingModel::hit_set_t hits;

	for( int binder = 0; 2 != binder; ++binder )
	{
		for( unsigned c = 0; 2 != c; ++c )
		{
			for( int p = 0; 2 != p; ++p )
			{
				for( unsigned l = 10; 30 != l; l += 10 )
				{
					for( double b = 0.0; 1.0 > b; b += 0.4 )
					{
						hits.insert(
							BindingModel::hit_t(
								0 == binder
									? BiobaseBindingModel::parameter_t( TableLink( MATRIX_DATA, 328 ) ).get_model()
									: BiobaseBindingModel::parameter_t( TableLink( MATRIX_DATA, 621 ) ).get_model(),
								b,
								p,
								l,
								1 == c ) );
					}
				}
			}
		}
	}

	check_serialisation( hits );
	check_strict_weak_ordering< BindingModel >( hits );
}


void
register_binding_model_tests(boost::unit_test::test_suite * test)
{
	test->add( BOOST_TEST_CASE( &check_ar_binding_model ), 0 );
	test->add( BOOST_TEST_CASE( &check_deaf_binding_model ), 0 );
	test->add( BOOST_TEST_CASE( &check_binding_hit_less_than ), 0 );
	test->add( BOOST_PARAM_TEST_CASE( &check_pssm_bayesian_binding_model, pssm_links.begin(), pssm_links.end() ), 0 );
	test->add( BOOST_PARAM_TEST_CASE( &check_biobase_binding_model, sequences.begin(), sequences.end() ), 0 );
}
