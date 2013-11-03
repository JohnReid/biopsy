header "pre_include_hpp" {
	#include "bio/factor.h"
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





class FactorParser extends BiobaseParser;
options {
	importVocab = BiobaseLexer;
}
{
	typedef Factor entry_t;
	
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
factor_name
	:
	{ sps.push_lexer("string"); }
	s:STRING { e->name = s->getText(); }
	{ sps.pop_lexer(); }
	;
	
protected
identifier : ignored_value ;

protected
binding_factor : ignored_value ;
	
protected
description : ignored_value ;
	
protected
gene 
	:
	{ sps.push_lexer("string"); }
	table_link_value[ e->gene ]
	{ sps.pop_lexer(); }
	ignored_value
	;
	
protected
sequence : ignored_value ;
	
protected
synonyms
	:
	{ sps.push_lexer("dot_separated"); }
	synonym_list_entry
	( ( SEMI | COLON ) synonym_list_entry )*
	DOT
	{ sps.pop_lexer(); }
	{ sps.push_lexer("string"); }
	( STRING )? //if a dot in the synonym list - disregard everything after it to the new line
	{ sps.pop_lexer(); }
	;

protected
synonym_list_entry
	:
	s:STRING { std::string res = s->getText(); } 
	(
		t:DOUBLECOLON { res += t->getText(); }
		u:STRING { res += u->getText(); }
	)*
	{e->synonyms.insert(res);}
	;

protected
taxonomy
	:
	{ sps.push_lexer("list"); }
	t1:STRING { e->taxonomies.push_back(t1->getText()); }
	(
		(
			SEMI
			t2:STRING { e->taxonomies.push_back(t2->getText()); }
		)+
		|
	)
	{ sps.pop_lexer(); }
	;


protected
matrices
{
	TableLink link;
}
	: 
	{ sps.push_lexer("list"); }
	table_link_value[link] { e->matrices.push_back( link ); }
	SEMI
	{ sps.pop_lexer(); }
	ignored_value
	;
	
protected
end_position //actually subunits for a factor not an end position at all
{
	TableLink link;
}
	: 
	{ sps.push_lexer("list"); }
	table_link_value[link] { e->subunits.push_back( link ); }
	SEMI
	{ sps.pop_lexer(); }
	ignored_value
	;
	
protected
complex
{
	TableLink link;
}
	: 
	{ sps.push_lexer("list"); }
	table_link_value[link] { e->complexes.push_back( link ); }
	SEMI
	{ sps.pop_lexer(); }
	ignored_value
	;
	
protected
super_family
{
	TableLink link;
}
	: 
	{ sps.push_lexer("list"); }
	table_link_value[link] { e->super_families.push_back( link ); }
	SEMI
	{ sps.pop_lexer(); }
	ignored_value
	;


protected
sub_family
{
	TableLink link;
}
	: 
	{ sps.push_lexer("list"); }
	table_link_value[link] { e->sub_families.push_back( link ); }
	{ sps.pop_lexer(); }
	ignored_value
	;

	
protected
type 
	: 
	{ sps.push_lexer("dot"); }
	s:STRING { e->_type = parse_factor_type( s->getText() ); }
	(DOT)?
	{ sps.pop_lexer(); }
	;

