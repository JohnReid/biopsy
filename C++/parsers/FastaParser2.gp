header "pre_include_hpp" {
	#include "bio/shared_parser_state.h"
	#include "bio/fasta.h"
}

options {
	language="Cpp";
	namespace = "bio";
}



class FastaParser2 extends Parser;
options {
	importVocab = FastaLexer2;
}
{
public:
	SharedParserState sps;
	fasta_map_t * fasta_map;
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