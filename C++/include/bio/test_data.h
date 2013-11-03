
#ifndef BIO_TEST_DATA
#define BIO_TEST_DATA

#include "bio/defs.h"
#include "bio/binding_model.h"
#include "bio/bifa_algorithm.h"
#include "bio/transcription_factor.h"

#include <boost/shared_ptr.hpp>

#include <vector>
#include <map>

BIO_NS_START


typedef double test_param_t;

struct ROCPoint
{
	double specificity;
	double sensitivity;

	ROCPoint( double specificity, double sensitivity );
};
std::ostream &
operator<<( std::ostream & os, const ROCPoint & roc_point );






struct BinaryTestResults
	: boost::addable< BinaryTestResults >
{
	unsigned true_positives;
	unsigned false_positives;
	unsigned true_negatives;
	unsigned false_negatives;

	BinaryTestResults(
		unsigned true_positives = 0,
		unsigned false_positives = 0,
		unsigned true_negatives = 0,
		unsigned false_negatives = 0);

	BinaryTestResults & operator()(bool tested_positive, bool should_be_positive);

	ROCPoint get_roc_point() const;

	BinaryTestResults & operator+=( const BinaryTestResults & rhs );
};





struct Test
{
	typedef boost::shared_ptr< Test > ptr_t;
	typedef std::vector< ptr_t > vec_t;

	virtual ~Test();

	virtual BinaryTestResults get_results( double threshold ) = 0;
};




BinaryTestResults
get_results_for(
	const BindingModel::hit_set_t & hits,
	double threshold,
	const BindingModel::set_t & universe,
	const BindingModel::set_t & true_positives);



struct BiFaTestData
{
public:
	typedef boost::shared_ptr< BiFaTestData > ptr_t;
	typedef std::vector< ptr_t > vec_t;
	typedef std::map< BiFaAlgorithm::ptr_t, BiFaOutput::ptr_t > output_map_t;
	typedef std::map< BiFaAlgorithm::ptr_t, TF::hit_set_ptr_t > factor_map_t;

protected:
	BiFaInput input;
	output_map_t output;
	factor_map_t factor_hits;

public:
	const BiFaInput & get_input() const;
	const TF::hit_set_t & get_factors_for( BiFaAlgorithm::ptr_t algorithm );

	BiFaTestData(
		const BiFaInput & bifa_input = BiFaInput() );

	virtual ~BiFaTestData();

	const BiFaOutput & get_output_for( BiFaAlgorithm::ptr_t algorithm );

	virtual Test::ptr_t get_test_for_algorithm( BiFaAlgorithm::ptr_t algorithm ) = 0;
};




struct BiFaFactorTest
	: Test
{
	BiFaTestData & data;
	const FactorLinkList & factors;
	BiFaAlgorithm::ptr_t algorithm;

	BiFaFactorTest(
		BiFaTestData & data,
		const FactorLinkList & factors,
		BiFaAlgorithm::ptr_t algorithm );

	virtual ~BiFaFactorTest();

	virtual BinaryTestResults get_results( double threshold );
};



typedef std::multimap< test_param_t, ROCPoint > test_result_map_t;
std::ostream &
operator<<( std::ostream & os, const test_result_map_t & result_map );



struct BiFaTestData2Test
{
	BiFaAlgorithm::ptr_t algorithm;

	BiFaTestData2Test( BiFaAlgorithm::ptr_t algorithm );

	Test::ptr_t operator()( BiFaTestData::ptr_t data ) const;
};



/**
Run tests.
*/
void
run_tests(
	Test::vec_t & test_data,
	test_result_map_t & results,
	double min_threshold,
	double max_threshold,
	unsigned num_thresholds );


void 
create_test_data_from_matrices(
	BiFaTestData::vec_t & test_data,
	size_t seq_length );

void 
create_test_data_from_fragments( 
	BiFaTestData::vec_t & test_data );

void 
create_test_data_from_sites( 
	BiFaTestData::vec_t & test_data );





BIO_NS_END

#endif //BIO_TEST_DATA
