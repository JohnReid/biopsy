/* Copyright John Reid 2007
*/

#include "bio-pch.h"

#include "bio/lexer.h"
#include "bio/shared_parser_state.h"

#include "BiobaseLexer.hpp"
#include "ArrayLexer.hpp"
#include "ColonStringLexer.hpp"
#include "ColonDotSemiStringLexer.hpp"
#include "CsvLexer.hpp"
#include "DotLexer.hpp"
#include "DotSemiStringLexer.hpp"
#include "DotSeparatedLexer.hpp"
#include "IdLexer.hpp"
#include "IndexLexer.hpp"
#include "LinkLexer.hpp"
#include "ListLexer.hpp"
#include "NameValuePairLexer.hpp"
#include "NumberLexer.hpp"
#include "PathwayLexer.hpp"
#include "SemiStringLexer.hpp"
#include "SequenceLexer.hpp"
#include "StringLexer.hpp"
#include "WhiteSpaceLexer.hpp"

#include <antlr/Parser.hpp>

BIO_NS_START

Lexer::Lexer( std::istream & is )
{
	add_lexer( "biobase",			lexer_ptr( new BiobaseLexer( is ) ) );

	add_lexer( "array",				lexer_ptr( new ArrayLexer( lexers[ "biobase" ]->getInputState() ) ) );
	add_lexer( "csv",				lexer_ptr( new CsvLexer( lexers[ "biobase" ]->getInputState() ) ) );
	add_lexer( "colondotsemi",		lexer_ptr( new ColonDotSemiStringLexer( lexers[ "biobase" ]->getInputState() ) ) );
	add_lexer( "colon",				lexer_ptr( new ColonStringLexer( lexers[ "biobase" ]->getInputState() ) ) );
	add_lexer( "dot",				lexer_ptr( new DotLexer( lexers[ "biobase" ]->getInputState() ) ) );
	add_lexer( "dot_separated",		lexer_ptr( new DotSeparatedLexer( lexers[ "biobase" ]->getInputState() ) ) );
	add_lexer( "dotsemi",			lexer_ptr( new DotSemiStringLexer( lexers[ "biobase" ]->getInputState() ) ) );
	add_lexer( "id",				lexer_ptr( new IdLexer( lexers[ "biobase" ]->getInputState() ) ) );
	add_lexer( "index",				lexer_ptr( new IndexLexer( lexers[ "biobase" ]->getInputState() ) ) );
	add_lexer( "list",				lexer_ptr( new ListLexer( lexers[ "biobase" ]->getInputState() ) ) );
	add_lexer( "link",				lexer_ptr( new LinkLexer( lexers[ "biobase" ]->getInputState() ) ) );
	add_lexer( "name_value_pair",	lexer_ptr( new NameValuePairLexer( lexers[ "biobase" ]->getInputState() ) ) );
	add_lexer( "number",			lexer_ptr( new NumberLexer( lexers[ "biobase" ]->getInputState() ) ) );
	add_lexer( "pathway",			lexer_ptr( new PathwayLexer( lexers[ "biobase" ]->getInputState() ) ) );
	add_lexer( "semi",				lexer_ptr( new SemiStringLexer( lexers[ "biobase" ]->getInputState() ) ) );
	add_lexer( "sequence",			lexer_ptr( new SequenceLexer( lexers[ "biobase" ]->getInputState() ) ) );
	add_lexer( "string",			lexer_ptr( new StringLexer( lexers[ "biobase" ]->getInputState() ) ) );
	add_lexer( "ws",				lexer_ptr( new WhiteSpaceLexer( lexers[ "biobase" ]->getInputState() ) ) );

	selector.select( "biobase" ); //start state
}

antlr::TokenStreamSelector & 
Lexer::get_selector() 
{ 
	return selector; 
}

void 
Lexer::add_lexer( const std::string & name, lexer_ptr l )
{
	if( lexers.end() != lexers.find( name ) )
	{
		throw std::logic_error( BIO_MAKE_STRING( "Already added lexer with name: " << name ) );
	}
	lexers[ name ] = l;
	selector.addInputStream( l.get(), name.c_str() );
}



void 
SharedParserState::push_lexer(const char * name) {
	if( antlr::DEBUG_PARSER ) std::cout << "Pushing lexer \"" << name << "\"\n";
	selector->push(name);
}
void 
SharedParserState::pop_lexer() {
	if( antlr::DEBUG_PARSER ) std::cout << "Popping lexer\n";
	selector->pop();
}
void 
SharedParserState::log(const char * text) const {
	std::cout << text << std::endl;
}	
void 
SharedParserState::log(const std::string & text) const {
	std::cout << text << std::endl;
}	
int 
SharedParserState::parse_number(const std::string & text) {
	return atoi(text.c_str());
}

float_t 
SharedParserState::parse_float(const std::string & text) {
	return (float_t)atof(text.c_str());
}


BIO_NS_END
