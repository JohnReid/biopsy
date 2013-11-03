header "pre_include_hpp" {
	#include "bio/site.h"
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





class SiteParser extends BiobaseParser;
options {
	importVocab = BiobaseLexer;
}
{
	typedef Site entry_t;
	
public:
	SharedParserState sps;
	entry_t::ptr_t e;
	entry_t::map_t * map;
	seq_t _seq; //store the sequence for the current site in here (sometimes this is split across 2 lines)
	
protected:
	void start_entry() {
		e.reset(new entry_t);
		_seq="";
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
			&& "" != e->sequence;
	}
}

protected
reference_point
	:
	{ sps.push_lexer("string"); }
	str_value[e->reference_point]
	{ sps.pop_lexer(); }
	;

protected
start_position
{ int value; }
	:
	value = integer_value
	{ e->start_position = value; }
	;

protected
end_position
{ int value; }
	:
	value = integer_value
	{ e->end_position = value; }
	;

