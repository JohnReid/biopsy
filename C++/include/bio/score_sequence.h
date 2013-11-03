
#ifndef BIO_SCORE_SEQUENCE_H_
#define BIO_SCORE_SEQUENCE_H_

#include "bio/defs.h"
#include "bio/match_hit.h"

BIO_NS_START



/** Runs the scorer over every position in the sequence and inserts the results into the result_insert_iterator. */
template <class Scorer, class SeqIt, class ResultInsIt>
void
score_sequence(
	const Scorer & scorer, //the scorer (e.g. pssm)
	SeqIt seq_begin, //the sequence start
	SeqIt seq_end, //the sequence end
	bool match_complement, //match the reversed complement?
	ResultInsIt result_insert_iterator) //put the match results here
{
	//for each sub-sequence
	for (int pos = 0; unsigned(seq_end - seq_begin) >= scorer.size(); ++seq_begin, ++pos)
	{
		*result_insert_iterator++ = typename ResultInsIt::container_type::value_type(scorer.score(seq_begin, false), pos, false);

		if (match_complement) {
			*result_insert_iterator++ = typename ResultInsIt::container_type::value_type(scorer.score(seq_begin, true), pos, true);
		}
	}
}


/** Runs the scorer over every position in the sequence and returns the best score. */
template <class Scorer, class SeqIt>
Hit
score_sequence(
	const Scorer & scorer,
	SeqIt seq_begin,
	SeqIt seq_end,
	bool match_complement)
{
	Hit result(-1.0, 0); //return -1.0 for error value if there is nothing to score over

	//for each sub-sequence
	typename Scorer::score_t score;
	for (int pos = 0; unsigned(seq_end - seq_begin) >= scorer.size(); ++seq_begin, ++pos) {
		score = scorer.score(seq_begin, false);
		if (score > result.score)
		{
			result.score = score;
			result.position = pos;
			result.complement = false;
		}

		if (match_complement)
		{
			score = scorer.score(seq_begin, true);
			if (score > result.score)
			{
				result.score = score;
				result.position = pos;
				result.complement = true;
			}
		}
	}

	return result;
}



BIO_NS_END


#endif //BIO_SCORE_SEQUENCE_H_
