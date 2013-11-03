header "pre_include_hpp" {
	#include "bio/matrix.h"	
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
}

options {
	language="Cpp";
	namespace = "bio";
}





class MatrixParser extends BiobaseParser;
options {
	importVocab = BiobaseLexer;
}
{
	typedef Matrix entry_t;
	
public:
	SharedParserState sps;
	entry_t::ptr_t e; //the matrix we are currently parsing
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
		return e->get_size() != 0 && UNKNOWN_DATA != e->accession_number.table_id;
	}
}


protected
name
	:
	{ sps.push_lexer("string"); }
	str_value[e->factor_name]
	{ sps.pop_lexer(); }
	;

protected
number_code
	: { sps.push_lexer("array"); }
	(SPACES)? na:STRING SPACES 
	nc:STRING SPACES 
	ng:STRING SPACES 
	nt:STRING SPACES 
	iupac:STRING
	{ sps.pop_lexer(); } 
	{ //could check here that the NUMBER_CODE is the correct index
		float_t counts[4];
		counts[0] = sps.parse_float(na->getText());
		counts[1] = sps.parse_float(nc->getText());
		counts[2] = sps.parse_float(ng->getText());
		counts[3] = sps.parse_float(nt->getText());
		e->pssm.push_back(PssmEntry(counts));

		e->consensus_matrix.push_back(iupac->getText()[0]);
	}
	;

protected
matrix_basis
	:
	{ sps.push_lexer("string"); }
	str_value[e->matrix_basis]
	{ sps.pop_lexer(); }
	{
		e -> number_of_sites = sps.parse_number(e->matrix_basis);
		if (e -> number_of_sites == 0)
		{
			if (strncmp(e->matrix_basis.c_str(),"computed",8) == 0)
				e -> number_of_sites = ANALOGUE_SITE_EQUIVALENT;
			else 
			{
				const char * p;
				if(p = strchr(e->matrix_basis.c_str(),':'))
				{
					double s = atof(p+1);
					if (s)
						e -> number_of_sites = (int)s;
				}
				if (e -> number_of_sites == 0)
				{
					if(strstr(e->matrix_basis.c_str(),"(SELEX)"))
						e -> number_of_sites = ANALOGUE_SITE_EQUIVALENT;
					else if(strstr(e->matrix_basis.c_str(),"Kd values"))
						e -> number_of_sites = ANALOGUE_SITE_EQUIVALENT;
					else if(strstr(e->matrix_basis.c_str(),"ChIP-chip"))
						e -> number_of_sites = ANALOGUE_SITE_EQUIVALENT;
					else
					{
						int a = 1;
					}
				}
			}
		}
	}
	;
protected
binding_site
{ AlignDescPtr ad(new AlignDesc); }
	: binding_site_value[*ad]
	{ e->align_descs.push_back(ad); }
	;
	
binding_site_value [ AlignDesc & ad ]
{ int matrix_start, matrix_length; }
	: { sps.push_lexer("list"); }		
	str_value[ad.sequence] SEMI //sequence segment
	{ sps.push_lexer("csv"); } table_link_value[ad.site] (COMMA table_link_value[ad.secondary_site])? { sps.pop_lexer(); } SEMI //binding site link number
	matrix_start = integer_value SEMI 
	matrix_length = integer_value SEMI 
	{ sps.push_lexer("csv"); } (
		g1:STRING { ad.gaps.clear(); ad.gaps.push_back(sps.parse_number(g1->getText())); }
		(COMMA g:STRING { ad.gaps.push_back(sps.parse_number(g->getText())); })*
	)? { sps.pop_lexer(); } SEMI 
	orientation:STRING
	{ sps.pop_lexer(); }
	{
		ad.start = matrix_start;
		ad.length = matrix_length;
		if ("p." == orientation->getText()) {
			ad.positive_orientation = true;
		} else if ("n." == orientation->getText()) {
			ad.positive_orientation = false;
		} else {
			throw std::logic_error( "Orientation must be '+'ve or '-'ve" );
		}
	}
	;	
	
protected
external_databases : ignored_value ;
	
protected
sequence : ignored_value ;
	

