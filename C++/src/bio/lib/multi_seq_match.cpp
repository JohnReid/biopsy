/* Copyright John Reid 2007
*/

#include "bio-pch.h"


#include "bio/defs.h"

#include "bio/multi_seq_match.h"
#include "bio/biobase_db.h"




BIO_NS_START

std::pair<bool, TableLink> get_common_factor(const TableLink & tl1, const TableLink & tl2)
{
	BiobaseTablePssmEntry * pssm_1 = BiobaseDb::singleton().get_pssm_entry(tl1);
	BiobaseTablePssmEntry * pssm_2 = BiobaseDb::singleton().get_pssm_entry(tl2);

	if (0 != pssm_1 && 0 != pssm_2)
	{
		for (FactorLinkList::const_iterator i1 = pssm_1->get_factors().begin();
			pssm_1->get_factors().end() != i1;
			++i1)
		{
			for (FactorLinkList::const_iterator i2 = pssm_2->get_factors().begin();
				pssm_2->get_factors().end() != i2;
				++i2)
			{
				//if the names match a case insensitive comparison that is good enough for us
				if (0 == STRICMP(i1->get()->name.c_str(), i2->get()->name.c_str()))
				{
					return std::make_pair(true, i1->get()->link);
				}
			}
		}
	}

	return std::make_pair(false, TableLink());
}

/** we say two hits are equal if they refer to a common factor. */
bool CommonFactorEqualTo::operator()(const TableLink & h1, const TableLink & h2) const
{
	return get_common_factor(h1, h2).first;
}


BIO_NS_END

