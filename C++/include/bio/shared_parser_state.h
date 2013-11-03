#ifndef BIOBASE_SHARED_PARSER_STATE_H_
#define BIOBASE_SHARED_PARSER_STATE_H_

#include "bio/defs.h"

#include <antlr/TokenStreamSelector.hpp>
#include <antlr/CharScanner.hpp>

#include <iostream>


BIO_NS_START

struct SharedParserState {
	antlr::TokenStreamSelector * selector;

/**	An attempt to manage the underlying stream.   The antlr people clearly intend
	this to work, but the documentation is poor, and there are no examples.  Also
	it is not easy to get at the underlying Input buffer to make it work
	
	'mark' where you want to go back to
	'rewind' to go back to it (which probaly means resetting the input state to 
		clear any current tokens)
	'commit' to throw the mark away.  The code looks as if it should keep the internal
		buffer down to a reasonable size after the commit, but it does not, so the
		extra reset stops uncontrolled buffer growth 
	*/

	unsigned int markInputstream() {return scanner().mark();};
	void rewindInputstream(unsigned int mark) {return scanner().rewind(mark);};
	void commitInputstream() {
		scanner().commit();
		scanner().getInputBuffer().reset();};


	void push_lexer( const char * name );
	void pop_lexer();
	void log( const char * text ) const;
	void log( const std::string & text ) const;
	static int parse_number( const std::string & text );
	static float_t parse_float( const std::string & text );
protected:
	antlr::CharScanner & scanner() 
		{return *static_cast<ANTLR_USE_NAMESPACE(antlr)CharScanner *> (selector -> getCurrentStream());};
};

BIO_NS_END


#endif //BIOBASE_SHARED_PARSER_STATE_H_
