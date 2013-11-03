header "pre_include_hpp"
{
	#include "bio/factor.h"
	#include "parser_includes.h"
}

header "pre_include_cpp"
{
	#include "bio-pch.h"
	#include "bio/defs.h"
}

header "post_include_cpp"
{
	using namespace std;
	using namespace antlr;
	using namespace BIO_NS;
	
#ifdef _MSC_VER
        #pragma warning( disable : 4101 )
#endif
}

options
{
	language="Cpp";
	namespace = "bio";
}





class DatabaseRefParser extends Parser;
options
{
	importVocab = Common;
	defaultErrorHandler = false;
}
{
public:
	Database db; /**< The database the reference refers to. */
	std::string acc; /**< The accession string. */
	bool parsed;
	std::set< std::string > unparsed_dbs;
	
	db_ref ref() const
	{
		return parse_db_ref_as( acc, db );
	}

protected:
	bool want_to_parse() const;
}

database_ref
	:
		database
		entry
	;

database
{ db = UNKNOWN_DB; }
	:
		d:STRING //{ sps.log( d->getText() ); } //match the colon delimited database
		{ db = parse_database_string( d->getText() ); }
		{ if( UNKNOWN_DB == db ) unparsed_dbs.insert( d->getText() ); }
		(COLON | DOT | SEMI)
	;
	
entry
{ parsed = false; }
	:
		{ want_to_parse() }? //do we want to parse this type of database reference?
		(
			//{ cout << db << ": want to parse\n"; }
			accession
			extra
			{ if( " _" != acc ) parsed = true; } //don't waste time with empty accessions
		) | (
			//we don't want to parse this entry
			//{ cout << db << ": don't want to parse\n"; }
			( STRING | DOT | SEMI | COLON )*
		)
	;
	
accession
{ acc = ""; }
	:
		( STRING DOT STRING )=> (
			s1:STRING DOT s2:STRING
			{ acc = BIO_MAKE_STRING( s1->getText() << "." << s2->getText() ); }
		) | (
			a:STRING 
			{ acc = a->getText(); }
		)
		{ acc = a->getText(); }
		//{ sps.log( acc->getText().c_str() ); } //match the semi delimited accession
	;
	
extra
	:
		( STRING | DOT | SEMI | COLON )+
	;

new_line
	:
		NEW_LINE
	;
