
#include "bio/contingency_homogeneity.h"
#include "bio/partition.h"
USING_BIO_NS

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/progress.hpp>
#include <boost/array.hpp>
#include <boost/assign/list_of.hpp>
using namespace boost;
using namespace boost::assign;
using boost::unit_test::test_suite;

#include <iostream>
using namespace std;


//#define VERBOSE


typedef ContingencyTable<2,2> ode_table_t;
ode_table_t ode_table;
ode_table_t::lambda_vec_t ode_lambda;

typedef ContingencyTable<2,25> mobility_table_t;
mobility_table_t mobility_table;
mobility_table_t::lambda_vec_t mobility_lambda;

void
build_data()
{
	static bool data_built = false;
	if (! data_built)
	{
		ode_table.data = 
			list_of
				(list_of( 9)(12))
				(list_of(78)(13));
		std::fill(ode_lambda.begin(), ode_lambda.end(), 1.0);

		mobility_table.data = 
			list_of
				(list_of( 50)( 18))
				(list_of( 45)( 17))
				(list_of(  8)( 16))
				(list_of( 18)(  4))
				(list_of(  8)(  2))
				(list_of( 28)( 24))
				(list_of(174)(105))
				(list_of( 84)(109))
				(list_of(154)( 59))
				(list_of( 55)( 21))
				(list_of( 11)( 23))
				(list_of( 78)( 84))
				(list_of(110)(289))
				(list_of(223)(217))
				(list_of( 96)( 95))
				(list_of( 14)(  8))
				(list_of(150)( 49))
				(list_of(185)(175))
				(list_of(714)(348))
				(list_of(447)(198))
				(list_of(  3)(  6))
				(list_of( 42)(  8))
				(list_of( 72)( 69))
				(list_of(320)(201))
				(list_of(411)(246))
				;
		std::fill(mobility_lambda.begin(), mobility_lambda.end(), 1.0);

		data_built = true;
	}
}

template <class PartitionIt>
void
print_partitions(PartitionIt it)
{
	typedef std::vector<char> data_set_t;
	static data_set_t v;

	// initialize the set
	if (v.empty())
	{
		for (unsigned i = 0; i < 26; ++i)
		{
			v.push_back('a' + i);
		}
	}

	try
	{
		while (true)
		{
			std::cout << it << " : " << it.subsets() << " : ";

			std::auto_ptr<std::vector<data_set_t> > part = it[v];
			std::cout << *part << '\n';

			++it;
		}
	}
	catch (std::overflow_error&) //how the partition code determines there are no more partitions!
	{
	}
}

template <class PartitionIt>
unsigned
get_num_partitions(PartitionIt it)
{
	unsigned result = 0;
	try
	{
		while (true)
		{
			++result;

			++it;
		}
	}
	catch (std::overflow_error&) //how the partition code determines there are no more partitions!
	{
	}
	return result;
}

void
check_partitions()
{
	cout << "******* check_partitions()" << endl;

	const unsigned size = 4;

#ifdef VERBOSE
	//all partitions
	print_partitions(partition::iterator(size));
	cout << endl << endl;

	//the partitions of size = 3
	print_partitions(partition::iterator_k(size, 3));
#endif //VERBOSE

	BOOST_CHECK_EQUAL(5u, get_num_partitions(partition::iterator(3)));
	BOOST_CHECK_EQUAL(1u, get_num_partitions(partition::iterator_k(3, 3)));
	BOOST_CHECK_EQUAL(3u, get_num_partitions(partition::iterator_k(3, 2)));
	BOOST_CHECK_EQUAL(1u, get_num_partitions(partition::iterator_k(3, 1)));

	BOOST_CHECK_EQUAL(15u, get_num_partitions(partition::iterator(4)));
	BOOST_CHECK_EQUAL(1u, get_num_partitions(partition::iterator_k(4, 4)));
	BOOST_CHECK_EQUAL(6u, get_num_partitions(partition::iterator_k(4, 3)));
	BOOST_CHECK_EQUAL(7u, get_num_partitions(partition::iterator_k(4, 2)));
	BOOST_CHECK_EQUAL(1u, get_num_partitions(partition::iterator_k(4, 1)));
}

void
check_contingency_calculate_L()
{
	cout << "******* check_calculate_L()" << endl;

	build_data();

	BOOST_CHECK_CLOSE(
		3.9831225868046277e-027, //not verified except by correct calculation of bayes factor for ode data
		contingency_calculate_L(
			1,
			ode_table.data.begin(),
			ode_table.data.end(),
			ode_lambda.begin(),
			ode_lambda.end()),
		0.001);
	BOOST_CHECK_CLOSE(
		2.3831634926329567e-024, //not verified except by correct calculation of bayes factor for ode data
		contingency_calculate_L(
			2,
			ode_table.data.begin(),
			ode_table.data.end(),
			ode_lambda.begin(),
			ode_lambda.end()),
		0.001);
}

void
print_contingency_calculate_gamma_exponential_prior(
	unsigned m,
	unsigned k,
	double mean)
{
	ExponentialPrior exp_prior(1.0);

	cout
		<< "m = " << m
		<< ", k = " << k
		<< ", mean = " << mean
		<< ", result = " << contingency_calculate_gamma(m, k, &exp_prior)
		<< endl;
}


void
check_contingency_calculate_gamma_exponential_prior()
{
	cout << "******* contingency_calculate_gamma()" << endl;

	build_data();

#ifdef VERBOSE
	print_contingency_calculate_gamma_exponential_prior(1, 2, 1);
	print_contingency_calculate_gamma_exponential_prior(2, 2, 1);
#endif //VERBOSE

	ExponentialPrior exp_prior(1.0);

	BOOST_CHECK_CLOSE(
		0.59634736,
		contingency_calculate_gamma(1, 2, &exp_prior),
		0.001);
	BOOST_CHECK_CLOSE(
		0.40365264,
		contingency_calculate_gamma(2, 2, &exp_prior),
		0.001);
	BOOST_CHECK_CLOSE(
		0.59634736,
		contingency_calculate_gamma(3, 2, &exp_prior),
		0.001);
	BOOST_CHECK_CLOSE(
		1.4036526,
		contingency_calculate_gamma(4, 2, &exp_prior),
		0.001);
}

template <class ContingencyTable>
void
t_check_contingency_calculate_bayes_factor(
	const ContingencyTable & table,
	const typename ContingencyTable::lambda_vec_t & lambda,
	const gamma_vec_t & contingency_gamma_exponential_prior,
	double expected_bayes_factor)
{
#ifdef VERBOSE
	boost::progress_timer timer;
#endif

	const double bayes_factor = 
		contingency_calculate_bayes_factor(
			contingency_gamma_exponential_prior.begin(),
			contingency_gamma_exponential_prior.end(),
			table.data.begin(),
			table.data.end(),
			lambda.begin(),
			lambda.end());

	const double bayes_factor_2 =
		contingency_calculate_L(1, table.data.begin(), table.data.end(), lambda.begin(), lambda.end())
			/ contingency_calculate_L(2, table.data.begin(), table.data.end(), lambda.begin(), lambda.end());

#ifdef VERBOSE
	cout << "Bayes factor: "
		<< bayes_factor
		<< endl;
	cout
		<< "L1(n) / L2(n): "
		<< bayes_factor_2
		<< endl;
#endif //VERBOSE

	BOOST_CHECK_CLOSE(bayes_factor, expected_bayes_factor, 0.001);
	BOOST_CHECK_CLOSE(bayes_factor_2, expected_bayes_factor, 0.001);
}

void
check_contingency_calculate_bayes_factor()
{
	cout << "******* check_contingency_calculate_bayes_factor()" << endl;

	build_data();

	ExponentialPrior exp_prior(1.0);

	gamma_vec_t contingency_gamma_exponential_prior;
	contingency_gamma_exponential_prior.push_back(
		contingency_calculate_gamma(1, 2, &exp_prior));
	contingency_gamma_exponential_prior.push_back(
		contingency_calculate_gamma(2, 2, &exp_prior));

	t_check_contingency_calculate_bayes_factor(
		ode_table,
		ode_lambda,
		contingency_gamma_exponential_prior,
		0.001671359);

	//takes too long and produces not a number
	if (false)
	{
		t_check_contingency_calculate_bayes_factor(
			mobility_table,
			mobility_lambda,
			contingency_gamma_exponential_prior,
			0.0);
	}
}

void register_contingency_homogeneity_tests(test_suite * test)
{
	test->add(BOOST_TEST_CASE(&check_contingency_calculate_gamma_exponential_prior), 0);
	test->add(BOOST_TEST_CASE(&check_partitions), 0);
	test->add(BOOST_TEST_CASE(&check_contingency_calculate_L), 0);
	test->add(BOOST_TEST_CASE(&check_contingency_calculate_bayes_factor), 0);
}


