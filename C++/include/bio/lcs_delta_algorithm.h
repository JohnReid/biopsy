
#ifndef BIO_LCS_DELTA_ALGORITHM_H_
#define BIO_LCS_DELTA_ALGORITHM_H_

#include "bio/defs.h"

#include <vector>
#include <functional>

BIO_NS_START

enum LcsResultElement
{
	LCS_ZERO,
	LCS_UP,
	LCS_LEFT,
	LCS_UP_LEFT,
	LCS_BANG
};

typedef std::vector<LcsResultElement> LcsResultVec;
typedef std::vector<LcsResultVec> LcsResultMat;

typedef std::vector<unsigned> LcsCountVec;
typedef std::vector<LcsCountVec> LcsCountMat;

/** Taken from http://en.wikipedia.org/wiki/Longest-common_subsequence_problem.
*/
template <class SeqIt, class EqualTo>
void
create_lcs_match_matrices(
	SeqIt x_begin,
	SeqIt x_end,
	SeqIt y_begin,
	SeqIt y_end,
	LcsResultMat & b,
	LcsCountMat & c,
	EqualTo equal_to = std::equal_to<typename SeqIt::value_type>())
{
	typedef typename SeqIt::difference_type diff_t;
	const diff_t m = x_end - x_begin;
	const diff_t n = y_end - y_begin;

	//resize the matrices
	b.resize(m + 1);
	c.resize(m + 1);
	for (diff_t i = 0; i != m + 1; ++i)
	{
		b[i].resize(n + 1);
		c[i].resize(n + 1);
	}

	//initialise
	for (diff_t i = 1; i != m + 1; ++i)
	{
		c[i][0] = 0;
		b[i][0] = LCS_ZERO;
	}
	for (diff_t j = 0; j != n + 1; ++j)
	{
		c[0][j] = 0;
		b[0][j] = LCS_ZERO;
	}

	SeqIt x = x_begin;
	for (diff_t i = 1; i != m + 1; ++i, ++x)
	{
		SeqIt y = y_begin;
		for (diff_t j = 1; j != n + 1; ++j, ++y)
		{
			if (equal_to(*x, *y))
			{
				c[i][j] = c[i-1][j-1]+1;
				b[i][j] = LCS_UP_LEFT;
			}
			else if (c[i-1][j] >= c[i][j-1])
			{
				c[i][j] = c[i-1][j];
				b[i][j] = LCS_UP;
			}
			else
			{
				c[i][j] = c[i][j-1];
				b[i][j] = LCS_LEFT;
			}
		}
    }
}

/**
Taken from http://en.wikipedia.org/wiki/Longest-common_subsequence_problem.
The current implementation inserts the reversed LCS into insert_it. To avoid this use rbegin and rend
instead of begin and end. 
This version of the algorithm allows us to use some sort of wildcarded equals operator and
retrieve the possibly different subsequences from the 2 input strings.
*/
template <class SeqIt, class Inserter, class EqualTo>
void
lcs_delta(
	SeqIt x_begin,
	SeqIt x_end,
	SeqIt y_begin,
	SeqIt y_end,
	Inserter inserter,
	EqualTo equal_to)
{
	typedef typename SeqIt::difference_type diff_t;

	LcsResultMat b;
	LcsCountMat c;

	create_lcs_match_matrices(
		x_begin,
		x_end,
		y_begin,
		y_end,
		b,
		c,
		equal_to);

	const diff_t m = b.size() - 1;
	const diff_t n = ((0 == m) ? 0 : b[0].size() - 1);

	diff_t i = m;
	diff_t j = n;
	while (0 != c[i][j])
	{
		switch (b[i][j])
		{
		case LCS_UP:
			--i;
			--x_end;
			break;

		case LCS_LEFT:
			--j;
			--y_end;
			break;

		case LCS_UP_LEFT:
			--i;
			--j;
			--x_end;
			--y_end;
			inserter(*x_end, *y_end);
			break;

		case LCS_ZERO:
		default:
			throw std::logic_error( "Should not get here" );
		}
	}
}

/**
Taken from http://en.wikipedia.org/wiki/Longest-common_subsequence_problem.
The current implementation inserts the reversed LCS into insert_it. To avoid this use rbegin and rend
instead of begin and end. 
*/
template <class SeqIt, class Inserter>
void
lcs_delta(
	SeqIt x_begin,
	SeqIt x_end,
	SeqIt y_begin,
	SeqIt y_end,
	Inserter inserter)
{
	lcs_delta(
		x_begin,
		x_end,
		y_begin,
		y_end,
		inserter,
		std::equal_to<typename SeqIt::value_type>());
}

/**
Stores both (potentially different due to wildcard equal operators) LCS sequences.
*/
template <class ContainerType>
struct LcsDoubleInserter
{
	typedef ContainerType container_t;
	typedef typename container_t::value_type value_t;

	container_t & x;
	container_t & y;

	LcsDoubleInserter(container_t & x, container_t & y) : x(x), y(y) { }

	void operator()(const value_t & cx, const value_t & cy)
	{
		x.push_back(cx);
		y.push_back(cy);
	}
};

/**
Stores just one (potentially different due to wildcard equal operators) LCS sequences.
*/
template <class ContainerType>
struct LcsSingleInserter
{
	typedef ContainerType container_t;
	typedef typename container_t::value_type value_t;

	container_t & x;

	LcsSingleInserter(container_t & x) : x(x) { }

	void operator()(const value_t & cx, const value_t & cy)
	{
		x.push_back(cx);
	}
};

template <class ContainerT>
LcsSingleInserter<ContainerT>
make_lcs_inserter(ContainerT & container)
{
	return LcsSingleInserter<ContainerT>(container);
}

template <class ContainerT>
LcsDoubleInserter<ContainerT>
make_lcs_inserter(ContainerT & container_1, ContainerT & container_2)
{
	return LcsDoubleInserter<ContainerT>(container_1, container_2);
}

BIO_NS_END



#endif //BIO_LCS_DELTA_ALGORITHM_H_
