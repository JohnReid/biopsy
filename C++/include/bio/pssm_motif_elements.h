#ifndef BIO_PSSM_MOTIF_ELEMENTS_H_
#define BIO_PSSM_MOTIF_ELEMENTS_H_

#include "bio/defs.h"
#include "bio/pssm_motif.h"
#include "bio/common.h"


BIO_NS_START



struct PssmMotifPssmElement : PssmMotif::ElementMatcher
{
	TableLink link;

	PssmMotifPssmElement(TableLink link);
	virtual ~PssmMotifPssmElement() { }

	bool matches(const MatchResults & match);
};


struct PssmMotifFactorElement : PssmMotif::ElementMatcher
{
	std::string factor_name;

	PssmMotifFactorElement(const std::string & factor_name);
	virtual ~PssmMotifFactorElement() { }

	bool matches(const MatchResults & match);
};


struct OrPssmMotifElement : PssmMotif::ElementMatcher
{
	PssmMotif::ElementMatcher::ptr_t e1;
	PssmMotif::ElementMatcher::ptr_t e2;

	OrPssmMotifElement(PssmMotif::ElementMatcher::ptr_t e1, PssmMotif::ElementMatcher::ptr_t e2);
	virtual ~OrPssmMotifElement() { }

	bool matches(const MatchResults & match);
};




BIO_NS_END



#endif //BIO_PSSM_MOTIF_ELEMENTS_H_
