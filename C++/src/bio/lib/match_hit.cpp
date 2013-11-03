/* Copyright John Reid 2007
*/

#include "bio-pch.h"


#include "bio/defs.h"


#include "bio/match_hit.h"

BIO_NS_START


/** Finds both scores (complementary and non-complementary) for a given position. Returns negative scores if not found. */
score_pair_t
find_scores_for_position(const hit_vec_t & hits, size_t position)
{
	score_pair_t result(-1.0, -1.0);

	hit_vec_t::const_iterator i = hits.begin();

	while (hits.end() != i)
	{
		if (size_t(i->position) == position)
		{
			if (i->complement)
			{
				result.second = i->score;
			}
			else
			{
				result.first = i->score;
			}
		}

		//did we find both scores?
		if (result.second >= 0.0 && result.first >= 0.0)
		{
			//yes
			break;
		}

		++i;
	}

	return result;
}

float_t
find_best_score_for_position(const hit_vec_t & hits, size_t position)
{
	const score_pair_t scores = find_scores_for_position(hits, position);
	return std::max(scores.first, scores.second);
}

size_t
rank_score(const hit_vec_t & hits, float_t score)
{
	size_t num_better = 0;
	size_t num_equal = 0;

	for (hit_vec_t::const_iterator i = hits.begin();
		hits.end() != i;
		++i)
	{
		if (score < i->score)
		{
			++num_better;
		}
		else if (score == i->score)
		{
			++num_equal;
		}
	}

	return num_better + num_equal / 2;
}

/** Output to stream. */
std::ostream &
operator<<(std::ostream & os, const Hit & hit)
{
	return os << hit.score << "," << hit.position << (hit.complement ? ",-" : ",+");
}

BIO_NS_END
