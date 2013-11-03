header "pre_include_hpp" {
	#include "bio/fragment.h"
	#include "parser_includes.h"
}

header "pre_include_cpp"
{
	#include "bio-pch.h"
	#include "bio/defs.h"
	#include "DatabaseRefParser.hpp"
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





class FragmentParser extends BiobaseParser;
options {
	importVocab = BiobaseLexer;
}
{
	typedef Fragment entry_t;
	
public:
	SharedParserState sps;
	entry_t::ptr_t e; //current entry
	entry_t::map_t * map;
	seq_t _seq; //store the sequence for the current site in here (sometimes this is split across 2 lines)
	
protected:
	void start_entry() {
		e.reset(new entry_t);
	}
	void end_entry() {
		if (is_entry_complete()) {
			(*map)[e->accession_number] = e;
		}
		e.reset();
	}
	bool is_entry_complete() {
		return
			UNKNOWN_DATA != e->accession_number.table_id;
	}
}

protected
identifier : ignored_value ;

protected
description // e.g. DE  Gene: G018345; Gene: G018461.
	:
		{ sps.push_lexer("colondotsemi"); }
		gene_reference ( SEMI gene_reference )* DOT
		{ sps.pop_lexer(); }
	;
	
protected
gene_reference
{ TableLink tl; }
	:
	STRING COLON table_link_value[tl]
	{ e->genes.push_back( tl ); }
	;
	
protected
external_databases : ignored_value ;


