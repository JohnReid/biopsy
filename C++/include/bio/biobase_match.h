#ifndef BIO_BIOBASE_MATCH_H_
#define BIO_BIOBASE_MATCH_H_

#include "bio/defs.h"
#include "bio/sequence.h"
#include "bio/match_hit.h"
#include "bio/matrix.h"
#include "bio/site.h"
#include "bio/cache.h"
#include "bio/singleton.h"

#include <algorithm>
#include <cmath>
#include <vector>

BIO_NS_START


/**
ConsT (consensus type) is often IupacCode or PssmEntry.
*/
template <class ConsT>
struct BiobasePssmEntry
	: public ConsT
{
	typedef ConsT consensus_t;

	BiobasePssmEntry(const consensus_t & e) : consensus_t(e)
	{
		//for each nucleotide
		conservation_information = 0.0;
		if (0 != this->get_num_observations())
		{
			for (std::string::const_iterator i = nucleotides.begin(); i != nucleotides.end(); ++i)
			{
				float_t fi = this->get_freq(*i);
				conservation_information += 0.0f == fi ? 0.0f : fi * std::log(4.0f * fi);
				if (BIO_ISNAN(conservation_information))
				{
					throw std::logic_error( "conservation information is NaN" );
				}
			}
		}
	}

	float_t conservation_information;
	float_t get_conservation_information() const { return conservation_information; }
};

/**
ConsT (consensus type) is often IupacCode or PssmEntry.
*/
template <class ConsT>
struct BiobasePssm
	: public std::vector<BiobasePssmEntry<ConsT> >
{
	typedef ConsT consensus_t;
	typedef BiobasePssmEntry<consensus_t> entry_t;
	typedef std::vector<entry_t> base_t;

	template <class ConsIt>
	BiobasePssm(ConsIt & begin, ConsIt & end) {
		max_score = 0.0;
		min_score = 0.0;
		for (ConsIt i = begin; i != end; ++i) {
			const entry_t e(*i);
			push_back(e);
			max_score += ((ConsT) *i).get_max() * e.conservation_information;
			min_score += ((ConsT) *i).get_min() * e.conservation_information;
		}
		if (BIO_ISNAN(max_score)) {
			throw std::logic_error( "max_score is NaN" );
		}
		if (BIO_ISNAN(min_score)) {
			throw std::logic_error( "min_score is NaN" );
		}
	}
	float_t max_score;
	float_t min_score;

	float_t get_range() const { return max_score - min_score; }
	float_t get_min_score() const { return min_score; }
};


Pssm make_pssm( TableLink link );
Pssm make_pssm( const Matrix * matrix );
Pssm make_pssm( Matrix::ptr_t matrix );
bool is_matchable(const Matrix * matrix);
bool is_matchable(const Matrix::ptr_t matrix);

template <class IupacIt>
Pssm
make_pssm_from_iupac(IupacIt seq_begin, IupacIt seq_end)
{
	Pssm result;
	for (IupacIt i = seq_begin;
		seq_end != i;
		++i)
	{
		IupacCode code(*i);

		float_t counts[4];
		counts[0] = (float_t)code.get_num('a');
		counts[1] = (float_t)code.get_num('c');
		counts[2] = (float_t)code.get_num('g');
		counts[3] = (float_t)code.get_num('t');

		result.push_back(PssmEntry(counts));
	}
	return result;
}

Pssm make_pssm( const Site * site );
Pssm make_pssm( const Site::ptr_t site );
bool is_matchable(const Site * site);
bool is_matchable(const Site::ptr_t site);



BIO_NS_END

#endif //BIO_BIOBASE_MATCH_H_
