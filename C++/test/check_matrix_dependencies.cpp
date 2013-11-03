
#include "bio/matrix_dependencies.h"
#include "bio/biobase_db.h"
#include "bio/biobase_filter.h"
#include "bio/biobase_data_traits.h"
USING_BIO_NS

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/progress.hpp>
#include <boost/assign/list_of.hpp>
using namespace boost;
using namespace boost::assign;
using boost::unit_test::test_suite;

#include <gsl/gsl_math.h>

#include <set>
#include <iostream>
using namespace std;


void
check_matrix_dependencies_arg(pssm_contingency_tables_t & contingency_tables)
{
	pssm_dependency_results_t results;
	calculate_dependencies(
		contingency_tables,
		results);

	//ordered set of results
	pssm_dependency_result_set result_set;
	build_result_set(results, result_set);

#ifdef VERBOSE
	unsigned max_to_show = 5;
	for (pssm_dependency_result_set::const_iterator i = result_set.begin();
		result_set.end() != i && 0 != max_to_show;
		++i, --max_to_show)
	{
		cout
			<< (*i)->first << " : " << (*i)->second.homogeneity_bayes_factor << '\n'
			<< contingency_tables[(*i)->first]
			<< '\n';
	}
#endif //VERBOSE
}

void
check_simple_contingency_dependencies_arg(base_pair_contingency_table_t first)
{
	//make the other entry in transpose of the first
	base_pair_contingency_table_t reverse = first;
	for (unsigned b1 = 0; 4 != b1; ++b1)
	{
		for (unsigned b2 = b1 + 1; 4 != b2; ++b2)
		{
			std::swap(reverse.data[b1][b2], reverse.data[b2][b1]);
		}
	}

	pssm_contingency_tables_t contingency_tables;
	contingency_tables[BasePair(0,1)] = first;
	contingency_tables[BasePair(1,0)] = reverse;

	check_matrix_dependencies_arg(contingency_tables);
}

void
check_simple_contingency_dependencies_table(base_pair_contingency_table_t::table_t data)
{
	base_pair_contingency_table_t first;
	first.data = data;
	check_simple_contingency_dependencies_arg(first);
}

void
check_simple_contingency_dependencies()
{
	cout << "******* check_simple_contingency_dependencies()" << endl;

	check_simple_contingency_dependencies_table(
		list_of
			(list_of(10)( 0)( 0)( 0))
			(list_of( 0)(10)( 0)( 0))
			(list_of( 0)( 0)(10)( 0))
			(list_of( 0)( 0)( 0)(10)));

	check_simple_contingency_dependencies_table(
		list_of
			(list_of(10)(10)(10)(10))
			(list_of(10)(10)(10)(10))
			(list_of(10)(10)(10)(10))
			(list_of(10)(10)(10)(10)));

	check_simple_contingency_dependencies_table(
		list_of
			(list_of(10)( 2)( 0)( 1))
			(list_of(10)( 2)( 0)( 1))
			(list_of(10)( 2)( 0)( 1))
			(list_of(10)( 2)( 0)( 1)));
}

struct MatrixResult
{
	const Matrix * matrix;
	BasePair base_pair;
	BasePairDependencyResults results;

	MatrixResult()
		: base_pair(0, 0)
	{
	}

	bool operator<(const MatrixResult & rhs) const
	{
		return results.homogeneity_bayes_factor < rhs.results.homogeneity_bayes_factor;
	}
};

struct MatrixResultLessThan
{
	double threshold;
	MatrixResultLessThan(double threshold) : threshold(threshold) { }
	bool operator()(const MatrixResult & result) { return result.results.homogeneity_bayes_factor < threshold; }
};

void
check_matrix_dependencies_arg(matrix_filter_it matrices_begin, matrix_filter_it matrices_end)
{
	//calculate the dependencies for all the matrices we can
	matrix_dependencies_map_t matrix_dependencies;
	std::set<MatrixResult> matrix_results;
	for (matrix_filter_it i = matrices_begin; 
		i != matrices_end;
		++i)
	{
		if (i->second->align_descs.empty())
		{
			continue;
		}

		try
		{
#ifdef VERBOSE_CHECKING
			cout << i->second->get_name() << '\n';
#endif

			matrix_dependencies[i->second.get()] = MatrixDependencies(i->second.get());
			
			for (pssm_dependency_results_t::const_iterator j = matrix_dependencies[i->second.get()].results.begin();
				matrix_dependencies[i->second.get()].results.end() != j;
				++j)
			{
				MatrixResult result;
				result.matrix = i->second.get();
				result.base_pair = j->first;
				result.results = j->second;
				matrix_results.insert(result);
			}
		}
#ifdef VERBOSE_CHECKING
		catch (const std::exception & ex)
		{
			cout << ex.what() << endl;
		}
		catch (const std::string & ex)
		{
			cout << ex.c_str() << endl;
		}
		catch (const char * ex)
		{
			cout << ex << endl;
		}
		catch (...)
		{
			cout << "Unknown exception" << endl;
		}
#else
		catch (...)
		{
		}
#endif
	}

#ifdef VERBOSE_CHECKING
	const double threshold = 0.09;
	cout
		<< "Gathered " << matrix_results.size() << " results; "
		<< std::count_if(matrix_results.begin(), matrix_results.end(), MatrixResultLessThan(threshold))
		<< " are less than " << threshold << "\n";
	unsigned max_to_show = 100;
	for (std::set<MatrixResult>::const_iterator i = matrix_results.begin();
		matrix_results.end() != i && 0 != max_to_show;
		++i, --max_to_show)
	{
		cout
			<< i->matrix->get_name() << '\t'
			<< i->base_pair << '\t'
			<< i->results.homogeneity_bayes_factor << '\n';
		cout << matrix_dependencies[i->matrix].contingency_tables[i->base_pair] << '\n';
	}
#endif //VERBOSE_CHECKING
}

void
check_matrix_dependencies()
{
	cout << "******* check_matrix_dependencies()" << endl;

	const matrix_filter_it matrices_begin(BiobaseDb::singleton().get_matrices().begin(), BiobaseDb::singleton().get_matrices().end());
	const matrix_filter_it matrices_end(BiobaseDb::singleton().get_matrices().end(), BiobaseDb::singleton().get_matrices().end());

	check_matrix_dependencies_arg(matrices_begin, matrices_end);
}

void
check_crebatf_dependencies()
{
	cout << "******* check_crebatf_dependencies()" << endl;

	Matrix::map_t::iterator begin = BiobaseDb::singleton().get_matrices().find(TableLink(MATRIX_DATA, 981));
	Matrix::map_t::iterator end = begin; ++end;

	const matrix_filter_it matrices_begin(begin, end);
	const matrix_filter_it matrices_end(end, end);

	check_matrix_dependencies_arg(matrices_begin, matrices_end);
}

void
check_simple_sequence_dependencies()
{
	cout << "******* check_simple_sequence_dependencies()" << endl;

	pssm_source_matrix_t pssm_source;
	add_sequence_to_pssm_source("ca", pssm_source);
	add_sequence_to_pssm_source("cc", pssm_source);
	add_sequence_to_pssm_source("cg", pssm_source);
	add_sequence_to_pssm_source("ct", pssm_source);

	pssm_contingency_tables_t contingency_tables;
	count_frequencies(
		pssm_source,
		contingency_tables);

	check_matrix_dependencies_arg(contingency_tables);
}



void register_matrix_dependencies_tests(test_suite * test)
{
	test->add(BOOST_TEST_CASE(&check_simple_sequence_dependencies), 0);
	test->add(BOOST_TEST_CASE(&check_simple_contingency_dependencies), 0);
	test->add(BOOST_TEST_CASE(&check_matrix_dependencies), 0);
	test->add(BOOST_TEST_CASE(&check_crebatf_dependencies), 0);
}


