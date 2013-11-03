/* Copyright John Reid 2007
*/

#include "bio-pch.h"


#include "bio/defs.h"
#include "bio/pssm.h"

#include <algorithm>
#include <numeric>



BIO_NS_START

PssmEntry::PssmEntry()
{
	for (size_t i = 0; i != 4; ++i)
	{
		counts[i] = 0;
		scores[i] = 0.0;
	}
}

PssmEntry::PssmEntry(
	float_t num_a,
	float_t num_c,
	float_t num_g,
	float_t num_t)
{
	counts[0] = num_a;
	counts[1] = num_c;
	counts[2] = num_g;
	counts[3] = num_t;

	init();
}

PssmEntry::PssmEntry(float_t * c)
{
	if( 0 == c ) {
		for (size_t i = 0; i != 4; ++i)
		{
			counts[i] = 0;
			scores[i] = 0.0;
		}
	} else {
		for (size_t i = 0; i != 4; ++i)
		{
			counts[i] = c[i];
		}
		init();
	}
}

float_t
PssmEntry::get_conservation_information() const
{
	//calculate the conservation information
	float_t conservation_information = 0.0;
	const float_t num_obs = get_num_observations();
	if (0 != num_obs)
	{
		//for each nucleotide
		for (size_t i = 0; i != 4; ++i)
		{
			float_t fi = ((float_t) counts[i]) / (float_t) num_obs;
			conservation_information +=
				(0.0f == fi)
					? 0.0f
					: fi * std::log(4.0f * fi);
			if (BIO_ISNAN(conservation_information))
			{
				throw std::logic_error( "conservation information is NaN" );
			}
		}
	}

	return conservation_information;
}

void
PssmEntry::init()
{
	const float_t conservation_information = get_conservation_information();

	//calculate the scores
	for (size_t i = 0; i != 4; ++i) //for each nucleotide
	{
		scores[i] = conservation_information * (float_t) counts[i];
	}
}


float_t
PssmEntry::get_num_observations() const
{
	return std::accumulate(counts, counts + 4, (float_t) 0);
}

float_t PssmEntry::get_max() const
{
	return *(std::max_element(scores, scores + 4));
}

float_t PssmEntry::get_min() const
{
	return *(std::min_element(scores, scores + 4));
}

float_t
PssmEntry::get_freq( char c, float_t pseudo_count ) const
{
	const float_t num_obs = float_t( get_num_observations() ) + float_t( 4 ) * pseudo_count;
	const float_t count = float_t( get_count(c) ) + pseudo_count;
	return
		( 0.0 == num_obs )
			? 0
			: count / num_obs;
}

Pssm::Pssm()
{
}


Pssm::~Pssm()
{
}


#if 0
bool pssm_print_debug_info = false;
#endif //0

float_t
Pssm::score(
	seq_t::const_iterator seq_begin,
	bool match_complement) const
{
#if 0
	if( pssm_print_debug_info )
	{
		std::cout << "Scoring pssm\n";
		std::cout << "this: " << this << "\n";
		std::cout << "Seq starts: " << *seq_begin << "\n";
	}
#endif //0

	float_t result = 0.0;
	float_t max = 0.0;
	float_t min = 0.0;

	// go backwards if matching complement
	if (match_complement)
	{
		for (const_reverse_iterator entry = rbegin();
			rend() != entry;
			++entry, ++seq_begin)
		{
			max += entry->get_max();
			min += entry->get_min();
			if ('n' == *seq_begin || 'N' == *seq_begin)
			{
				result += entry->get_min();
			}
			else
			{
				result += entry->get_score(complement(*seq_begin));
			}
		}
	}
	else
	{
		for (const_iterator entry = begin();
			end() != entry;
			++entry, ++seq_begin)
		{
			max += entry->get_max();
			min += entry->get_min();
			if ('n' == *seq_begin || 'N' == *seq_begin)
			{
				result += entry->get_min();
			}
			else
			{
				result += entry->get_score(*seq_begin);
			}
		}
	}

	//normalise into [0,1]
	const float_t range = max - min;
	result = (float_t)((0.0 == range) ? 0.0 : (result - min) / range);
	result = std::max(0.0f, result);
	result = std::min(result, 1.0f);

#if 0
	if( pssm_print_debug_info )
	{
		std::cout << "Scored pssm\n";
	}
#endif //0

	return result;
}


BIO_NS_END
