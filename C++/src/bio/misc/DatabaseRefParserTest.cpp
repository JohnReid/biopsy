/* Copyright John Reid 2007
*/

#include "bio-pch.h"



#include "DatabaseRefParser.hpp"
#include "ColonDotSemiStringLexer.hpp"
#include "bio/lexer.h"
#include "bio/counter.h"
#include <antlr/ANTLRException.hpp>
#include <antlr/Parser.hpp>

using namespace std;
using namespace antlr;
using namespace boost;
USING_BIO_NS;

void
parse_molecule_dr_stream( istream & in_stream )
{
	string line;
	while( getline( in_stream, line ) )
	{
		parse_molecule_dr_line( line );
	}
}

void
parse_molecule_dr_file()
{
	const string filename = "c:\\data\\biobase\\transfac\\molecule_dr.dat";
	cout << "Parsing file \"" << filename << "\"\n";
	ifstream is( filename.c_str() );

	parse_molecule_dr_stream( is );
}

int
main( int argc, char * argv[] )
{
	int result = 0;

	//antlr::DEBUG_PARSER = true;

	//string filename = "c:\\data\\biobase\\transfac\\dr_small.dat";
	string filename = "c:\\data\\biobase\\transfac\\dr.dat";
	cout << "Parsing file \"" << filename << "\"\n";
	ifstream is( filename.c_str() );

	ColonDotSemiStringLexer lexer( is );
	DatabaseRefParser parser( lexer );

	Counter< Database > counter;

	try
	{
		parse_molecule_dr_file();

		while( true )
		{
			parser.database_ref();
			if( parser.parsed )
			{
				const db_ref ref = parser.ref();
				counter.increment( ref.db );
				//cout << "Parsed: " << ref << "\n";
				if( UNKNOWN_DB == ref.db ) throw std::logic_error( "Badly parsed reference" );
			}
			else
			{
				//cout << "Not parsed: " << parser.db << "\n";
			}
			parser.new_line();
		}
	}
	catch (const ANTLRException & ex)
	{
		cerr
			<< to_simple_string( boost::posix_time::second_clock::local_time() )
			<< ": ANTLR Error: "
			<< ex.toString()
			<< endl;
		result = -1;
	}
	catch (const exception & ex)
	{
		cerr 
			<< to_simple_string( boost::posix_time::second_clock::local_time() )
			<< ": Error: " 
			<< ex.what() 
			<< endl;
		result = -1;
	}
	catch (const string & msg)
	{
		cerr 
			<< to_simple_string( boost::posix_time::second_clock::local_time() )
			<< ": Error: " 
			<< msg 
			<< endl;
		result = -2;
	}
	catch (const char * msg)
	{
		cerr 
			<< to_simple_string( boost::posix_time::second_clock::local_time() )
			<< ": Error: " 
			<< msg 
			<< endl;
		result = -3;
	}
	catch (...)
	{
		cerr 
			<< to_simple_string( boost::posix_time::second_clock::local_time() )
			<< ": Undefined error" 
			<< endl;
		result = -4;
	}

	counter.print();

	BOOST_FOREACH( const std::string & db, parser.unparsed_dbs ) std::cout << "Unparsed: " << db << "\n";

	return result;
}
