header "pre_include_hpp" {
	#include "bio/evidence.h"
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





class EvidenceParser extends BiobaseParser;
options {
	importVocab = BiobaseLexer;
}
{
public:
	typedef Evidence entry_t;

	SharedParserState sps;
	entry_t::ptr_t e;
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
		return UNKNOWN_DATA != e->accession_number.table_id;
	}
}

protected
identifier : ignored_value ;

protected
binding_factor : ignored_value ;
	
protected
composite_element
	:
	{ sps.push_lexer("link"); } 
	table_link_value[e->composite_element]
	{ sps.pop_lexer(); }
	ignored_value
	;

protected
description : ignored_value ;
	
protected
sequence : ignored_value ;
	

