/* Copyright John Reid 2007
*/

#include "bio-pch.h"


#include "bio/defs.h"

#include "bio/matrix_dependencies.h"
#include "bio/biobase_db.h"
#include "bio/biobase_filter.h"
#include "bio/biobase_data_traits.h"


BIO_NS_START

BasePair::BasePair(unsigned observed_index, unsigned conditioned_index)
	: observed_index(observed_index)
	, conditioned_index(conditioned_index)
{
}

bool
BasePair::operator<(const BasePair & rhs) const
{
	return
		observed_index < rhs.observed_index
		|| (observed_index == rhs.observed_index
			&& conditioned_index < rhs.conditioned_index);
}

std::ostream &
operator<<(std::ostream & os, const BasePair & base_pair)
{
	return os << base_pair.observed_index << '|' << base_pair.conditioned_index;
}

void
calculate_dependencies(
	pssm_contingency_tables_t & contingency_tables,
	pssm_dependency_results_t & results)
{
	using namespace boost::assign;

	ExponentialPrior exp_prior(1.0);

	//calculate the gamma vec once
	static const gamma_vec_t contingency_gamma_exponential_prior =
		list_of
 			(contingency_calculate_gamma(1, 4, &exp_prior))
			(contingency_calculate_gamma(2, 4, &exp_prior))
			(contingency_calculate_gamma(3, 4, &exp_prior))
			(contingency_calculate_gamma(4, 4, &exp_prior))
			;

	//initialise the lambda vec once
	static const base_pair_contingency_table_t::lambda_vec_t lambda = list_of(1.0)(1.0)(1.0)(1.0);

	//generate results for each contingency table
	results.clear();
	for (pssm_contingency_tables_t::iterator t = contingency_tables.begin();
		contingency_tables.end() != t;
		++t)
	{
		BasePairDependencyResults r;
		r.homogeneity_bayes_factor = 
			contingency_calculate_bayes_factor(
				contingency_gamma_exponential_prior.begin(),
				contingency_gamma_exponential_prior.end(),
				t->second.data.begin(),
				t->second.data.end(),
				lambda.begin(),
				lambda.end());
		results[t->first] = r;
	}
}


inline
unsigned
get_index_for(char c)
{
	switch(c)
	{
	case 'a': case 'A': return 0;
	case 'c': case 'C': return 1;
	case 'g': case 'G': return 2;
	case 't': case 'T': return 3;
	default:
		throw std::logic_error( "Unknown nucleotide" );
	}
}
	
void
count_frequencies(
	const pssm_source_matrix_t & pssm_source,
	pssm_contingency_tables_t & contingency_tables)
{
	using namespace std;

	const unsigned seq_length = pssm_source.size();
	const unsigned num_seqs = (0 != seq_length) ? pssm_source[0].size() : 0;

	//calculate the conditional distributions for pairs of positions
	//for each pair of positions, i and j
	for (unsigned i = 0; seq_length != i; ++i)
	{
		for (unsigned j = 0; seq_length != j; ++j)
		{
			//don't compare a position with itself
			if (j == i)
			{
				continue;
			}

			base_pair_contingency_table_t & table = contingency_tables[BasePair(i, j)];
			for (unsigned b1 = 0; 4 != b1; ++b1)
			{
				for (unsigned b2 = 0; 4 != b2; ++b2)
				{
					table.data[b1][b2] = 0;
				}
			}
			
			for (unsigned s = 0; num_seqs != s; ++s)
			{
				++table.data[get_index_for(pssm_source[i][s])][get_index_for(pssm_source[j][s])];
			}
		}
	}
}

void
build_result_set(
	pssm_dependency_results_t & results,
	pssm_dependency_result_set & result_set)
{
	result_set.clear();
	for (pssm_dependency_results_t::const_iterator i = results.begin();
		results.end() != i;
		++i)
	{
		result_set.insert(i);
	}
}

bool
HomogeneityLessThen::operator()(pssm_dependency_result_it lhs, pssm_dependency_result_it rhs) const
{
	if (lhs->second.homogeneity_bayes_factor < rhs->second.homogeneity_bayes_factor)
	{
		return true;
	}
	else if (lhs->second.homogeneity_bayes_factor == rhs->second.homogeneity_bayes_factor)
	{
		return lhs->first < rhs->first;
	}
	else
	{
		return false;
	}
}


template <class SeqIt>
pssm_source_matrix_t & 
add_seq_to_pssm_source(
	SeqIt seq_begin,
	SeqIt seq_end,
	pssm_source_matrix_t & pssm_source)
{
	typedef SeqIt seq_it;

	unsigned seq_length = 0;
	for (SeqIt s = seq_begin;
		seq_end != s;
		++s)
	{
		++seq_length;
	}

	//if we are the first sequence to be added...
	if (pssm_source.empty() && 0 != seq_length)
	{
		pssm_source.resize(seq_length);
	}
	else if (! pssm_source.empty() && pssm_source.size() != seq_length)
	{
		throw std::logic_error( "Sequences of unequal length" );
	}

	for (pssm_source_matrix_t::iterator p = pssm_source.begin();
		seq_end != seq_begin;
		++p, ++seq_begin)
	{
		if (pssm_source.end() == p)
		{
			throw std::logic_error( "Different length sequences" );
		}

		p->push_back(*seq_begin);
	}

	return pssm_source;
}

void
build_pssm_source(
	const Matrix * matrix,
	pssm_source_matrix_t & pssm_source)
{
	pssm_source.clear();

	for (AlignDescList::const_iterator i = matrix->align_descs.begin();
		matrix->align_descs.end() != i;
		++i)
	{
		//check it is not an artificial sequence...
		const Site * site = BiobaseDb::singleton().get_entry<SITE_DATA>((*i)->site);
		if (0 != site && "AS" != site->id.species_group)
		{
			add_seq_to_pssm_source((*i)->sequence.begin(), (*i)->sequence.end(), pssm_source);
		}
	}
}

void
add_sequence_to_pssm_source(
	const seq_t & sequence,
	pssm_source_matrix_t & pssm_source)
{
	add_seq_to_pssm_source(sequence.begin(), sequence.end(), pssm_source);
}


MatrixDependencies::MatrixDependencies()
{
}

MatrixDependencies::MatrixDependencies(const pssm_source_matrix_t & pssm_source)
: pssm_source(pssm_source)
{
	count_frequencies(pssm_source, contingency_tables);
	calculate_dependencies(contingency_tables, results);
	build_result_set(results, result_set);
}

MatrixDependencies::MatrixDependencies(const Matrix * matrix)
{
	build_pssm_source(matrix, pssm_source);
	count_frequencies(pssm_source, contingency_tables);
	calculate_dependencies(contingency_tables, results);
	build_result_set(results, result_set);
}

/** Designed to be functor in std::transform over Biobase pssm map. */
struct PssmSourceBuilder
{
	template< typename T >
	pssm_source_map_t::value_type operator()(const T & value)
	{
		pssm_source_map_t::value_type result = std::make_pair( value.first, pssm_source_matrix_ptr_t( new pssm_source_matrix_t ) );

		try
		{
			build_pssm_source( value.second.get(), *(result.second) );
		}
		catch (...)
		{
			std::cout << "Could not build pssm source for " << value.first << "\n";
			result.second.reset();
		}

		return result;
	}
};

struct PssmSourceIsEmpty
{
	bool operator()(pssm_source_map_t::value_type value) const
	{
		return 0 == value.second || 0 == value.second->size();
	}
};

void
build_all_pssm_sources(pssm_source_map_t & pssm_source_map)
{
	pssm_source_map_t tmp;

	std::transform(
		get_matrices_begin(),
		get_matrices_end(),
		std::inserter(tmp, tmp.begin()),
		PssmSourceBuilder());

	std::remove_copy_if(
		tmp.begin(),
		tmp.end(),
		std::inserter(pssm_source_map, pssm_source_map.begin()),
		PssmSourceIsEmpty());
}


std::ostream &
operator<<(std::ostream & os, pssm_source_map_t::value_type value)
{
	os << value.first << "\t" << BiobaseDb::singleton().get_pssm_entry(value.first)->get_name() << "\n";

	if (0 == value.second)
	{
		os << "No source matrix\n";
	}
	else
	{
		const pssm_source_matrix_t & source = *(value.second);
		for (unsigned seq = 0; ! source.empty() && source[0].size() != seq; ++seq)
		{
			for (unsigned pos = 0; source.size() != pos; ++pos)
			{
				os << source[pos][seq];
			}
			os << "\n";
		}
	}

	return os;
}


BIO_NS_END
