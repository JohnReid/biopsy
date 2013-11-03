/* Copyright John Reid 2007
*/

#include "bio-pch.h"


#include "bio/defs.h"

#include "bio/pssm_motif_elements.h"
#include "bio/biobase_data_traits.h"

BIO_NS_START




PssmMotifPssmElement::PssmMotifPssmElement(TableLink link)
	: link(link)
{
}

bool
PssmMotifPssmElement::matches(const MatchResults & match)
{
	return match.link == link;
}

PssmMotifFactorElement::PssmMotifFactorElement(const std::string & factor_name)
	: factor_name(bio_to_upper(factor_name))
{
	//check that the name matches at least one factor in biobase...

	for (Factor::map_t::const_iterator f = BiobaseDb::singleton().get_factors().begin();
		BiobaseDb::singleton().get_factors().end() != f;
		++f)
	{
		if (f->second->is_synonym(factor_name))
		{
			return;
		}
	}

	throw BIO_MAKE_STRING("Could not find any factor in Biobase that matches \"" << factor_name << "\"");
}

bool
PssmMotifFactorElement::matches(const MatchResults & match)
{
	//check if each factor is a synonym

	const FactorLinkList & factors = BiobaseDb::singleton().get_pssm_entry(match.link)->get_factors();
	for (FactorLinkList::const_iterator f = factors.begin();
		factors.end() != f;
		++f)
	{
		Factor * factor = BiobaseDb::singleton().get_entry<FACTOR_DATA>(f->get()->link);
		if (0 != factor && factor->is_synonym(factor_name))
		{
			return true;
		}
	}

	return false;
}

OrPssmMotifElement::OrPssmMotifElement(PssmMotif::ElementMatcher::ptr_t e1, PssmMotif::ElementMatcher::ptr_t e2)
	: e1(e1)
	, e2(e2)
{
}

bool
OrPssmMotifElement::matches(const MatchResults & match)
{
	return e1->matches(match) || e2->matches(match);
}


BIO_NS_END

