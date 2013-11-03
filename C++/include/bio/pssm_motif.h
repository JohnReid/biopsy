#ifndef BIO_PSSM_MOTIF_H_
#define BIO_PSSM_MOTIF_H_

#include "bio/defs.h"
#include "bio/run_match.h"

#include <boost/shared_ptr.hpp>

#include <string>
#include <vector>


BIO_NS_START



struct Score
{
	typedef double value_t;

protected:
	value_t score;

public:
	Score()
	{
		reset();
	}

	void reset()
	{
		score = 1.0;
	}

	void add(value_t evidence)
	{
		score *= (1.0 - evidence);
	}

	value_t get() const
	{
		return 1.0 - score;
	}
};





/** A PssmMotif can be thought of as a regular expression for searching in 
BiFa analyses. */
struct PssmMotif
{
	typedef boost::shared_ptr< PssmMotif > ptr_t;
	typedef std::vector< ptr_t > vec_t;
	typedef std::set< ptr_t > set_t;
	typedef std::map< ptr_t, Score > score_map_t;

	/** A range of distances permissible between consecutive elements of the motif. */
	struct Distance
	{
		typedef boost::shared_ptr< Distance > ptr_t;

		int min;
		int max;

		Distance(int min = -1, int max = -1)
			: min(min)
			, max(max)
		{
		}
	};

	/** One element of the PssmMotif. Has a virtual boolean method that decides whether a match
	occurs. */
	struct ElementMatcher
	{
		typedef boost::shared_ptr< ElementMatcher > ptr_t;

		virtual ~ElementMatcher() { }

		/** Does the Pssm match this motif element? */
		virtual bool matches(const MatchResults & match) = 0;
	};

	typedef std::pair< Distance::ptr_t, ElementMatcher::ptr_t > Element;
	typedef std::vector< Element > ElementVec;

	/** Describes where an Element matched in the BiFa analysis. */
	struct HitElement
	{
		typedef std::vector< HitElement > vec_t;

		match_result_vec_t::const_iterator match_result;
		ElementMatcher::ptr_t element;

		HitElement(
			match_result_vec_t::const_iterator match_result,
			ElementMatcher::ptr_t element)
			: match_result(match_result)
			, element(element)
		{
		}

		int get_end() const;
	};
	typedef HitElement::vec_t Hit;
	typedef std::vector< Hit > HitVec;
	static float_t score(const Hit & hit);

	ElementVec elements;

	static ptr_t parse(const std::string & description);

	static double get_score(const score_map_t & score_map);

	void find_in(const match_result_vec_t & matches, HitVec & hits);
};



std::ostream & operator<<(std::ostream & os, const PssmMotif::HitElement::vec_t & hit);

BIO_NS_END



#endif //BIO_PSSM_MOTIF_H_
