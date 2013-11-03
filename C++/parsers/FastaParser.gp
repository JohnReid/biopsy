header "pre_include_hpp" {
	#include "bio/shared_parser_state.h"
	#include <sstream>
}

header "pre_include_cpp"
{
	#include "bio-pch.h"
	#include "bio/defs.h"
}

options {
	language="Cpp";
	namespace = "bio";
}



class FastaParser extends Parser;
options {
	importVocab = FastaLexer;
}
{
public:
	SharedParserState sps;
	std::stringstream * sequence;
}

fasta
	:
	DESC_BEGIN DESC_LINE
	{ sps.push_lexer("nucleo"); }
	(
		n:NUCLEO_SEQ { *sequence << n->getText(); }
		NEW_LINE
	)+
	;