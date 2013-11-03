/* Copyright John Reid 2007
*/

#include "bio-pch.h"


#include "bio/defs.h"

#include "bio/pssm_motif.h"
#include "bio/biobase_db.h"

#include "PssmMotifLexer.hpp"
#include "PssmMotifParser.hpp"

BIO_NS_START

void
find_pssm_motif_hits(
	match_result_vec_t::const_iterator current_match,
	match_result_vec_t::const_iterator end_match,
	PssmMotif::ElementVec::const_iterator next_element_to_match,
	PssmMotif::ElementVec::const_iterator end_element,
	PssmMotif::Hit hit_so_far,
	PssmMotif::HitVec & hits)
{
	//did we match all the elements?
	if (end_element == next_element_to_match)
	{
		//yes so add the hit so far to the hits
		hits.push_back(hit_so_far);

		//all done
		return;
	}

	//for each match
	for (match_result_vec_t::const_iterator m = current_match;
		end_match != m;
		++m)
	{
		//do we already have one element matched? and is this hit within the range?
		if (! hit_so_far.empty())
		{
			//yes - we need to check the position...


			const PssmMotif::Distance::ptr_t distance = next_element_to_match->first;
			if (0 != distance)
			{
				//is it less than the minimum?
				if (m->result.position - hit_so_far.rbegin()->get_end() < distance->min)
				{
					//ignore and carry on
					continue;
				}

				//is it more than the maximum?
				if (m->result.position - hit_so_far.rbegin()->get_end() > distance->max)
				{
					//won't find any more in the correct range
					break;
				}
			}
		}

		//does it match the element?
		if (next_element_to_match->second->matches(*m))
		{
			//yes

			//add to the hit so far
			hit_so_far.push_back(PssmMotif::HitElement(m, next_element_to_match->second));

			//recurse - carry on finding the next element in the next matches
			find_pssm_motif_hits(
				m + 1,
				end_match,
				next_element_to_match + 1,
				end_element,
				hit_so_far,
				hits);

			hit_so_far.erase(hit_so_far.end() - 1);
		}
	}
}

int PssmMotif::HitElement::get_end() const
{
	return match_result->result.position + BiobaseDb::singleton().get_pssm_entry(match_result->link)->get_size();
}

void PssmMotif::find_in(const match_result_vec_t & matches, HitVec & hits)
{
	find_pssm_motif_hits(
		matches.begin(),
		matches.end(),
		elements.begin(),
		elements.end(),
		Hit(),
		hits);
}

double
PssmMotif::get_score(const score_map_t & score_map)
{
	double result = 1.0;

	for (score_map_t::const_iterator m = score_map.begin();
		score_map.end() != m;
		++m)
	{
		result *= m->second.get();
	}

	return result;
}

PssmMotif::ptr_t PssmMotif::parse(const std::string & description)
{
	std::stringstream stream(description);
	PssmMotifLexer lexer(stream);
	PssmMotifParser parser(lexer);

	lexer.found_eof = false;
	PssmMotif::ptr_t result = parser.pssm_motif();
	if (! lexer.found_eof)
	{
		throw std::logic_error( "Did not consume all input" );
	}

	return result;
}

float_t
PssmMotif::score(const PssmMotif::Hit & hit)
{
	float_t score = float_t(1.0);
	for (PssmMotif::Hit::const_iterator e = hit.begin();
		hit.end() != e;
		++e)
	{
		score *= (float_t(1.0) - e->match_result->result.score);
	}
	return float_t(1.0) - score;
}

std::ostream & operator<<(std::ostream & os, const PssmMotif::HitElement::vec_t & hit)
{
	for (PssmMotif::HitElement::vec_t::const_iterator e = hit.begin();
		hit.end() != e;
		++e)
	{
		os << '\t' << *(e->match_result) << '\n';
	}
	os << '\n';

	return os;
}


BIO_NS_END



