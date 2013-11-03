header "pre_include_hpp" {
	#include "bio/gene.h"
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





class GeneParser extends BiobaseParser;
options {
	importVocab = BiobaseLexer;
}
{
	typedef Gene entry_t;
	
public:
	SharedParserState sps;
	entry_t::ptr_t e; //current entry
	entry_t::map_t * map;
	
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
identifier
	: { sps.push_lexer("id"); }
	str_value[ e->species ]
	DOLLAR 
	( 
		str_value[ e->name ]
		(
			e1:UNDER{ e->name.append(e1->getText()); }
			| e2:STRING{ e->name.append(e2->getText()); }
		)*
	)?
	{ sps.pop_lexer(); }
	;

protected
description : ignored_value ;

protected
binding_factor : ignored_value ;

protected
sequence : ignored_value ;

protected
short_description : ignored_value ;

