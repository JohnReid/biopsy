header "pre_include_hpp" {
	#include "bio/matrix_match.h"	

	#include <fstream>	
}

header "pre_include_cpp"
{
	#include "bio-pch.h"
	#include "bio/defs.h"
}

header "post_include_hpp" {
	#include "MatrixMatchLexer.hpp"
}

options {
	language="Cpp";
	namespace = "bio";
}




class MatrixMatchLexer extends Lexer;
options {
	charVocabulary = '\3'..'\377';
	testLiterals=true;
}
tokens {
	AC="AC";
	M="M";
	ID="ID";
	NA="NA";
	MATR_LENGTH="MATR_LENGTH";
	CORE_START="CORE_START";
	CORE_LENGTH="CORE_LENGTH";
	MAXIMAL="MAXIMAL";
	MINIMAL="MINIMAL";
	THRESHOLD="THRESHOLD";
	WEIGHTS="WEIGHTS";
	A="A";
	C="C";
	G="G";
	T="T";
	SLASHES="//";
	NUMBER_OF_MATRICES="NUMBER_OF_MATRICES";
}
protected ALPHA			: ('A'..'Z') | ('a'..'z') ;
protected DIGIT			: ('0'..'9') ;
protected MINUS			: '-' ;
COLON					: ':' ;
STRING					: (ALPHA | '(' | '/') (~('\r' | '\n' | ' ' | ':'))* ; //starts with a char then matches everything up to whitespace or colon
NUMBER					: (MINUS)? (DIGIT | '.')+ ;
WS					    : (
							' '
							| '\t'
							| ( ('\r')? '\n' {newline();} )
						)+ { $setType(antlr::Token::SKIP); } ; //skip whitespace

class MatrixMatchParser extends Parser;
options {
}
{
protected:
	typedef MatrixMatch entry_t;
	typedef entry_t::map_t map_t;
	
	map_t * map;
	entry_t e; //the matrix we are currently parsing
	
	static int parse_int(const std::string & number)
	{
		return atoi(number.c_str());
	}
	static float_t parse_float(const std::string & number)
	{
		return (float_t) atof(number.c_str());
	}
	
public:
	static void parse(const std::string & filename, map_t & m)
	{
		std::ifstream is(filename.c_str());
		if (! is) {

			throw std::string("Could not open file: ") + filename;
		
		} else {
		
			MatrixMatchLexer lexer(is);
			MatrixMatchParser parser(lexer);
			
			parser.map = &m;

			parser.table();
		}
	}
}


table
	:
	(
		matrix_match { map->operator[](e.accession_number) = e; }
		SLASHES
	)+
	NUMBER_OF_MATRICES NUMBER
	;
	
protected
matrix_match
	:
	AC ac:STRING { e.accession_number = TableLink(MATRIX_DATA, parse_int(ac->getText().substr(1))); }
	ID STRING
	NA (~(MATR_LENGTH))+
	MATR_LENGTH NUMBER
	CORE_START NUMBER
	CORE_LENGTH NUMBER
	MAXIMAL NUMBER
	MINIMAL NUMBER
	THRESHOLD th:NUMBER { e.threshold = parse_float(th->getText()); }
	WEIGHTS
	(
		NUMBER
		A COLON NUMBER
		C COLON NUMBER
		G COLON NUMBER
		T COLON NUMBER
	)+
	;
	