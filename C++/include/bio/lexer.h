#ifndef BIOBASE_LEXER_H_
#define BIOBASE_LEXER_H_

#include "bio/defs.h"

#include <antlr/TokenStreamSelector.hpp>
#include <antlr/CharScanner.hpp>

BIO_NS_START


class Lexer 
{
public:
	Lexer( std::istream & is );

	antlr::TokenStreamSelector & get_selector();

protected:
	typedef boost::shared_ptr< ANTLR_USE_NAMESPACE(antlr)CharScanner > lexer_ptr;
	typedef std::map< std::string, lexer_ptr > lexer_map;

	lexer_map lexers;
	antlr::TokenStreamSelector selector;

	void add_lexer( const std::string & name, lexer_ptr l );
};

#if 0
class Lexer {
public:
	Lexer(std::istream & is)
		: biobaseLexer(is)
		, arrayLexer(biobaseLexer.getInputState())
		, csvLexer(biobaseLexer.getInputState())
		, dotLexer(biobaseLexer.getInputState())
		, dotSeparatedLexer(biobaseLexer.getInputState())
		, idLexer(biobaseLexer.getInputState())
		, indexLexer(biobaseLexer.getInputState())
		, listLexer(biobaseLexer.getInputState())
		, linkLexer(biobaseLexer.getInputState())
		, nameValuePairLexer(biobaseLexer.getInputState())
		, numberLexer(biobaseLexer.getInputState())
		, pathwayLexer(biobaseLexer.getInputState())
		, sequenceLexer(biobaseLexer.getInputState())
		, stringLexer(biobaseLexer.getInputState())
	{
		selector.addInputStream(&biobaseLexer, "biobase");
		selector.select("biobase"); //start state

		selector.addInputStream(&arrayLexer, "array");
		selector.addInputStream(&csvLexer, "colon");
		selector.addInputStream(&csvLexer, "csv");
		selector.addInputStream(&dotSeparatedLexer, "dot_separated");
		selector.addInputStream(&dotLexer, "dot");
		selector.addInputStream(&idLexer, "id");
		selector.addInputStream(&indexLexer, "index");
		selector.addInputStream(&listLexer, "list");
		selector.addInputStream(&linkLexer, "link");
		selector.addInputStream(&nameValuePairLexer, "name_value_pair");
		selector.addInputStream(&numberLexer, "number");
		selector.addInputStream(&pathwayLexer, "pathway");
		selector.addInputStream(&sequenceLexer, "sequence");
		selector.addInputStream(&stringLexer, "string");
	}

	antlr::TokenStreamSelector & get_selector() { return selector; }

protected:
	BiobaseLexer					biobaseLexer;
	ArrayLexer						arrayLexer;
	CsvLexer						csvLexer;
	DotLexer						dotLexer;
	DotSeparatedLexer				dotSeparatedLexer;
	IdLexer							idLexer;
	IndexLexer						indexLexer;
	ListLexer						listLexer;
	LinkLexer						linkLexer;
	NameValuePairLexer				nameValuePairLexer;
	NumberLexer						numberLexer;
	PathwayLexer					pathwayLexer;
	SequenceLexer					sequenceLexer;
	StringLexer						stringLexer;
	antlr::TokenStreamSelector		selector;
};
#endif



BIO_NS_END


#endif //BIOBASE_LEXER_H_
