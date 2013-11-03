header "pre_include_hpp" {
	#include "bio/compel.h"
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





class CompelParser extends BiobaseParser;
options {
	importVocab = Common;
}
{
public:
	typedef Compel entry_t;
	
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
binding_factor : ignored_value ;
	
protected
description : ignored_value ;
	
protected
binding_site
{ CompelBindingSite bs; }
	: binding_site_value[bs]
	{ e->binding_sites.push_back(bs); }
	;

/**	A potentially useful routine to skip to the next line when a field is split across
	multiple labelled lines */
skip_newline [int type] :
	{sps.push_lexer("biobase");
	match(type);
	match(SPACES); 
	sps.pop_lexer();} ;
	exception
	catch [std::logic_error const& ex] { sps.pop_lexer(); }
	
	
binding_site_value [ CompelBindingSite & bs ]
{ 	std::string position; }	:
    { sps.push_lexer("list"); }		
	str_value[position] SEMI //start to end
	str_value[bs.factor] SEMI //factor name(s) as string
	{ sps.pop_lexer(); } //list

	{ sps.push_lexer("link"); } //	Need to make sure that a tablelink is a valid alternative
	( NEW_LINE {skip_newline(BS);} )?
	{ sps.pop_lexer(); } //list
	
	table_link_value[bs.site_link] //link to compel site
	ignored_value
	{
		static const boost::regex position_regex( " *(-?[0-9]+) +to +(-?[0-9]+)" );
		boost::cmatch what;
		if( ! regex_match( position.c_str(), what, position_regex ) ) {
			throw std::logic_error( BIO_MAKE_STRING( "Could not parse \""<<position<<"\"" ) );
		}
		bs.start = detail::safe_lexical_cast< int >( what[1] );
		bs.end = detail::safe_lexical_cast< int >( what[2] );
	}
	;
	
protected
evidence
{ TableLink ev; }
	:
	{ sps.push_lexer("link"); } 
	table_link_value[ev] { e->evidences.push_back(ev); }
	{ sps.pop_lexer(); } 
	;
	
protected
gene_link
	:
	table_link_value[ e->gene ]
	;
	exception
	catch [std::logic_error const& ex] { sps.pop_lexer(); }

/** Ignore second line of multiline gene descriptions */	
protected
gene 
	:
	{ if (e->gene.entry_idx == 0) 
		{sps.push_lexer("string");
		gene_link();
		sps.pop_lexer();}
	}
	ignored_value
	;

protected
sequence //C00427 has more than one line of sequence
{ std::string s; }
	:
	{ sps.push_lexer("string"); }
	str_value[s]
	{ sps.pop_lexer(); }
	{ e->sequence += s; }
	;

protected
type
{ std::string s; }
	:
	{ sps.push_lexer("string"); }
	str_value[ s ]
	{ sps.pop_lexer(); }
	{
		if( s == "synergism" ) e->type = COMPEL_SYNERGISM;
		else if( s == "antagonism" ) e->type = COMPEL_ANTAGONISM;
		else throw std::logic_error( BIO_MAKE_STRING( "Cannot parse \""<<s<<"\" as compel type." ) );
	}
	;

protected
comment
{ std::string line; }
	: 
	{ sps.push_lexer("string"); }
	str_value[ line ]
	{ sps.pop_lexer(); }
	{ 
		e->comment += " "; 
		e->comment += line; 
	}
	;

protected
position 
{ std::string position; }
	:
	{ sps.push_lexer("string"); }
	str_value[position] //start to end
	{ sps.pop_lexer(); }
	{
		static const boost::regex position_regex( " *(-?[0-9]+) +to +(-?[0-9]+)" );
		boost::cmatch what;
		if( ! regex_match( position.c_str(), what, position_regex ) ) {
			throw std::logic_error( BIO_MAKE_STRING( "Could not parse \""<<position<<"\"" ) );
		}
		e->begin = detail::safe_lexical_cast< int >( what[1] );
		e->end = detail::safe_lexical_cast< int >( what[2] );
	}
	;

