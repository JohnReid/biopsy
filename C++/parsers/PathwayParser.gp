header "pre_include_hpp" {
	#include "bio/pathway.h"
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
}

options {
	language="Cpp";
	namespace = "bio";
}





class PathwayParser extends BiobaseParser;
options {
	importVocab = BiobaseLexer;
}
{
	typedef Pathway entry_t;
	
public:
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
		return
			UNKNOWN_DATA != e->accession_number.table_id
			&& UNKNOWN_PW != e->pathway_type;
	}
}

protected
identifier : ignored_value ;

protected
binding_factor : ignored_value ;
	
protected
description : ignored_value ;
	
protected
sequence : ignored_value ;
	
protected
external_databases : ignored_value ;
	
protected
type
	:
	{ sps.push_lexer("string"); }
	s:STRING
	{
		if (s->getText() == "chain.") {
			e->pathway_type = CHAIN_PW;
		} else if (s->getText() == "evidence chain.") {
			e->pathway_type = EVIDENCE_CHAIN_PW;
		} else if (s->getText() == "pathway.") {
			e->pathway_type = PATHWAY_PW;
		} else {
			e->pathway_type = UNKNOWN_PW;
		}
	}
	{ sps.pop_lexer(); }
	;

protected
super_family
	:
	super_family_value[e->super_families]
	;
	
protected //name of the pathway
name
	:
	{ sps.push_lexer("string"); }
	s:STRING { e->name = s->getText(); }
	{ sps.pop_lexer(); }
	;
	

