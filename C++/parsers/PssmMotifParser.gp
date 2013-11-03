header "pre_include_hpp" {
	#include "bio/pssm_motif.h"
	#include "bio/pssm_motif_elements.h"
}

header "pre_include_cpp"
{
	#include "bio-pch.h"
	#include "bio/defs.h"
}

header "post_include_cpp" {
	using namespace std;
	using namespace antlr;
	using namespace BIO_NS;
	
#ifdef _MSC_VER
	#pragma warning( disable : 4101 )
#endif
}

options {
	language="Cpp";
	namespace = "bio";
}





class PssmMotifParser extends Parser;
options {
	importVocab = PssmMotifLexer;
	defaultErrorHandler = false;
}

pssm_motif
	returns [ PssmMotif::ptr_t result ]
{
	result.reset(new PssmMotif);
	PssmMotif::ElementMatcher::ptr_t e;
	PssmMotif::Distance::ptr_t d;
}
	:
	e = motif_expression
	{ result->elements.push_back(std::make_pair(PssmMotif::Distance::ptr_t(), e)); }
	(
		SEMI
		( WS )?
		d = distance
		e = motif_expression
		{ result->elements.push_back(std::make_pair(d, e)); }
	)*
	;
	
motif_expression
	returns [ PssmMotif::ElementMatcher::ptr_t result ]
{
	PssmMotif::ElementMatcher::ptr_t e;
}
	:
	result = motif_element
	(
		OR
		WS
		e = motif_element
		{ result.reset(new OrPssmMotifElement(result, e)); }
	)*
	;

motif_element
	returns [ PssmMotif::ElementMatcher::ptr_t result ]
	:
	(
		result = pssm_motif_element
		|
		result = factor_motif_element
	)
	( WS )?
	;
	
pssm_motif_element
	returns [ PssmMotif::ElementMatcher::ptr_t result ]
{
	TableLink link;
}
	:
	PSSM
	WS
	link = biobase_link
	{
		result.reset(new PssmMotifPssmElement(link));
	}
	;
	
factor_motif_element
	returns [ PssmMotif::ElementMatcher::ptr_t result ]
	:
	FACTOR
	WS
	s:STRING
	{
		result.reset(new PssmMotifFactorElement(s->getText().substr(1, s->getText().length() - 2)));
	}
	;
	
distance
	returns [ PssmMotif::Distance::ptr_t result ]
{
	int min, max;
}
	:
	(
		OPEN_RANGE
		min = integer
		COMMA
		max = integer
		CLOSE_RANGE
		{
			result.reset(new PssmMotif::Distance(min, max));
		}
		( WS )?
	)?
	;
	
integer
	returns [ int result ]
	:
	i:INTEGER
	{
		result = atoi(i->getText().c_str());
	}
	;
	
biobase_link
	returns [ TableLink result ]
	:
	l:BIOBASE_LINK
	{
		result = parse_table_link_accession_number(l->getText());
	}
	;
	
