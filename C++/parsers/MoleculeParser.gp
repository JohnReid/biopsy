header "pre_include_hpp" {
	#include "bio/molecule.h"
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





class MoleculeParser extends BiobaseParser;
options {
	importVocab = BiobaseLexer;
}
{
	typedef Molecule entry_t;
	
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
			UNKNOWN_DATA != e->accession_number.table_id;
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
external_databases_unneeded
	:
		{ sps.push_lexer( "string" ); }
		s:STRING
		{ sps.pop_lexer(); }
		{ 
			const std::pair< bool, db_ref > parsed = parse_transfac_db_ref_line( s->getText() );
			if( parsed.first )
			{
				e->database_refs.push_back( parsed.second );
			}
		}
	;

//PW  <CH000000666>; alpha IIb beta3 ---> Rac1 (chain).
protected
pathway
{ TableLink l; }
	:
	{ sps.push_lexer("pathway"); }
	LESS
	{ sps.push_lexer("link"); }
	table_link_value[l] { e->pathways.push_back(l); }
	{ sps.pop_lexer(); }
	GREATER
	{ sps.pop_lexer(); }
	{ sps.push_lexer("string"); }
	STRING
	{ sps.pop_lexer(); }
	;
	
protected
super_family
	:
	super_family_value[e->super_families]
	;
	

