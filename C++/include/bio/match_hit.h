#ifndef BIO_MATCH_HIT_H_
#define BIO_MATCH_HIT_H_

#include "bio/defs.h"

#include <vector>
#include <list>


BIO_NS_START

/** Stores a score and a position where the match was. */
struct Hit
{
	Hit(
		float_t score = -1.0,
		int position = -1,
		bool complement = false)
		: score(score)
		, position(position)
		, complement(complement)
	{ }

	float_t score;
	int position;
	bool complement;

	bool operator==(const Hit & rhs) const
	{
		return position == rhs.position && complement == rhs.complement && score == rhs.score;
	}

	bool operator<(const Hit & rhs) const
	{
		if (position == rhs.position)
		{
			if (complement == rhs.complement)
			{
				return score < rhs.score;
			}
			else
			{
				return complement < rhs.complement;
			}
		}
		else
		{
			return position < rhs.position;
		}
	}

protected:
	friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & score;
		ar & position;
		ar & complement;
    }
};

/** An array of Hits. */
typedef std::list<Hit> hit_vec_t;

/** Output to stream. */
std::ostream &
operator<<(std::ostream & os, const Hit & hit);

/** A pair of scores from matching the sequence and its complement. */
typedef std::pair<float_t, float_t> score_pair_t;

/** Finds both scores (complementary and non-complementary) for a given position. Returns negative scores if not found. */
score_pair_t
find_scores_for_position(const hit_vec_t & hits, size_t position);

/** Finds best score (complementary or non-complementary) for a given position. Returns negative scores if not found. */
float_t
find_best_score_for_position(const hit_vec_t & hits, size_t position);

/** Ranks the score. I.e. what index it would be in an ordered list of the hits. Ties are averaged. I.e. a tie for
positions 3-7 would give a rank of 5. */
size_t
rank_score(const hit_vec_t & hits, float_t score);


BIO_NS_END

#endif //BIO_MATCH_HIT_H_
