
#ifndef BIO_MULTI_SEQ_MATCH_H_
#define BIO_MULTI_SEQ_MATCH_H_

#include "bio/defs.h"
#include "bio/run_match.h"
#include "bio/lcs_delta_algorithm.h"

BIO_NS_START

/** Get the common factor two matrices or sites represent. The bool result is false if there is none. */
std::pair<bool, TableLink> get_common_factor(const TableLink & tl1, const TableLink & tl2);

/** we say two hits are equal if they refer to a common factor. */
struct CommonFactorEqualTo
{
	bool operator()(const TableLink & h1, const TableLink & h2) const;
};

struct FactorInserter
{
	FactorInserter(TableLinkVec & factors)
		: factors(factors)
	{
	}

	void operator()(const TableLink & h1, const TableLink & h2)
	{
		std::pair<bool, TableLink> result = get_common_factor(h1, h2);
		if (result.first)
		{
			//std::cout << "Found match: " << result.second << std::endl;
			factors.push_back(result.second);
		}
	};

	TableLinkVec & factors;
};

template <class SeqsIt>
void
multiple_sequence_match(
	SeqsIt seqs_begin,
	SeqsIt seqs_end,
	float_t threshold,
	ScoreAlgorithm algorithm)
{
	typedef std::map<SeqsIt, TableLinkVec> hits_map_t;
	hits_map_t hits_map;

	//score each sequence
	BiobasePssmFilter filter;
	for (SeqsIt s1 = seqs_begin;
		seqs_end != s1;
		++s1)
	{

		match_result_vec_t hits;
		pssm_match(
			s1->first,
			s1->second,
			threshold,
			algorithm,
			false,
			filter,
			filter,
			std::inserter(hits, hits.begin()));

		//these are sorted by link rather than position - rearrange
		std::sort(hits.begin(), hits.end(), MatchResultsPositionLessThan());
		
		for (match_result_vec_t::const_iterator i = hits.begin();
			hits.end() != i;
			++i)
		{
			hits_map[s1].push_back(i->link);
		}
	}

	//for each pair of sequences
	for (SeqsIt s1 = seqs_begin;
		seqs_end != s1;
		++s1)
	{
		for (SeqsIt s2 = s1 + 1;
			seqs_end != s2;
			++s2)
		{
			//match the factors the hits represent
			TableLinkVec matched_factors;
			lcs_delta(
				hits_map[s1].rbegin(),
				hits_map[s1].rend(),
				hits_map[s2].rbegin(),
				hits_map[s2].rend(),
				FactorInserter(matched_factors),
				CommonFactorEqualTo());

			//match the hits
			TableLinkVec matched_hits;
			lcs_delta(
				hits_map[s1].rbegin(),
				hits_map[s1].rend(),
				hits_map[s2].rbegin(),
				hits_map[s2].rend(),
				make_lcs_inserter(matched_hits));

#if 0
			std::cout
				<< "Comparing "
				<< hits_map[s1].size()
				<< " hits against "
				<< hits_map[s2].size()
				<< " found "
				<< matched_factors.size()
				<< " matched factors and "
				<< matched_hits.size()
				<< " matched hits"
				<< std::endl;
#endif
		}
	}

	//the common longest common subsequence we find
	SeqsIt s1 = seqs_begin;
	if (seqs_end != s1)
	{
		TableLinkVec multiple_lcs = hits_map[s1];
		++s1;
		while (seqs_end != s1)
		{
			TableLinkVec result;
			lcs_delta(
				hits_map[s1].rbegin(),
				hits_map[s1].rend(),
				multiple_lcs.rbegin(),
				multiple_lcs.rend(),
				make_lcs_inserter(result));
			multiple_lcs = result;
			++s1;
		}
#if 0
		std::cout << "Found a longest common subsequence of length " << multiple_lcs.size() << std::endl;
		std::copy(multiple_lcs.begin(), multiple_lcs.end(), std::ostream_iterator<TableLink>(std::cout, ", "));
		std::cout << std::endl;
#endif
	}
}

BIO_NS_END


#endif //BIO_MULTI_SEQ_MATCH_H_
