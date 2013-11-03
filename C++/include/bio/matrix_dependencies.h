
#ifndef BIO_MATRIX_DEPENDENCIES_H_
#define BIO_MATRIX_DEPENDENCIES_H_

#include "bio/defs.h"
#include "bio/matrix.h"
#include "bio/contingency_homogeneity.h"

#include <boost/array.hpp>

#include <gsl/gsl_math.h>

#include <vector>
#include <numeric>
#include <set>

BIO_NS_START

/** Indexes a pair of bases for comparison. */
struct BasePair
{
	/** The index of the base we wish to analyse. */
	unsigned observed_index;

	/** The index of the base we condition on. */
	unsigned conditioned_index;

	BasePair(unsigned observed_index, unsigned conditioned_index);

	bool operator<(const BasePair & rhs) const;
};

std::ostream &
operator<<(std::ostream & os, const BasePair & base_pair);

/** A contingency table to compare a pair of bases. */
typedef ContingencyTable<4, 4> base_pair_contingency_table_t;

/** The contingency tables for all base pairs in a pssm. */
typedef std::map<BasePair, base_pair_contingency_table_t> pssm_contingency_tables_t;

/** The results of analysing the dependencies between a pair of bases. */
struct BasePairDependencyResults
{
	double homogeneity_bayes_factor;
};

/** The results of the dependencies between all base pairs in a pssm. */
typedef std::map< BasePair, BasePairDependencyResults > pssm_dependency_results_t;

/** An iterator pointing to a result. */
typedef pssm_dependency_results_t::const_iterator pssm_dependency_result_it;

/** A less than operator that lets us sort by scores. */
struct HomogeneityLessThen
{
	bool operator()(pssm_dependency_result_it lhs, pssm_dependency_result_it rhs) const;
};

/** An ordered set of results. */
typedef std::set< pssm_dependency_result_it, HomogeneityLessThen > pssm_dependency_result_set;

/** Contains all the sequences that were used to compose a Pssm. Indexed first by position then by sequence number. */
typedef std::vector< std::vector<char> > pssm_source_matrix_t;
typedef boost::shared_ptr< pssm_source_matrix_t > pssm_source_matrix_ptr_t;
typedef std::map< TableLink, pssm_source_matrix_ptr_t > pssm_source_map_t;


std::ostream &
operator<<(std::ostream & os, pssm_source_map_t::value_type value);

void
build_all_pssm_sources(pssm_source_map_t & pssm_source_map);

/** Take a Transfac matrix and generate a matrix from the sequences that defined it. */
void
build_pssm_source(
	const Matrix * matrix,
	pssm_source_matrix_t & pssm_source);

/** Adds a sequence to a pssm source. */
void
add_sequence_to_pssm_source(
	const seq_t & sequence,
	pssm_source_matrix_t & pssm_source);

/** Count the conditional frequencies in a matrix and generate contingency tables. */
void
count_frequencies(
	const pssm_source_matrix_t & pssm_source,
	pssm_contingency_tables_t & contingency_tables);

/** Calculate the dependencies in all the contingency tables. */
void
calculate_dependencies(
	pssm_contingency_tables_t & contingency_tables,
	pssm_dependency_results_t & results);

/** Build an ordered set of results. */
void
build_result_set(
	pssm_dependency_results_t & results,
	pssm_dependency_result_set & result_set);



struct MatrixDependencies
{
	pssm_source_matrix_t pssm_source;
	pssm_contingency_tables_t contingency_tables;
	pssm_dependency_results_t results;
	pssm_dependency_result_set result_set;

	MatrixDependencies(const pssm_source_matrix_t & pssm_source);
	MatrixDependencies(const Matrix * matrix);
	MatrixDependencies();
};

typedef std::map<const Matrix *, MatrixDependencies> matrix_dependencies_map_t;

BIO_NS_END




#endif //BIO_MATRIX_DEPENDENCIES_H_
