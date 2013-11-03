options {
	language="Cpp";
	namespace = "bio";
}





class BiobaseParser extends Parser;
options {
	importVocab = BiobaseLexer;
	k = 2;
}

table
	: (options {greedy=true;} : entry)+
	;
	
entry
	: { start_entry(); }
	(options {greedy=true;} : group)+
	END_OF_ENTRY NEW_LINE
	{ }
	{ end_entry(); }
	;
	
group
	: (options {greedy=true;} : value NEW_LINE)+
	(END_OF_GROUP NEW_LINE)?
	;

/** Switch between all possible field values. */
value
{ }
	: (
		CODE SPACES
		ignored_value
	) | (
		AC SPACES //Accession number
		accession_number
	) | (
		BA SPACES 	//Basis of data, often the number of samples used
		matrix_basis
	) | (
		BF SPACES 	//Binding factors 	  	factors shown to bind to this site (linked accession number; name; "Quality" of the factor-site interaction on a six level scale; biological species of the factor.)
		binding_factor
	) | (
		BS SPACES 	//Binding sites  	   	list of aligned sequence segments used for matrix generation (if available) followed by a link to the respective binding site in TRANSFAC (site accession number) or TRANSCompel from which the segment was derived and a description how it was derived (start, length, gaps and orientation of the depicted sequence segment - flanking gaps included - relative to the sequence in the site entry)
		binding_site
	) | (
		CC SPACES 	//comment?
		comment
	) | (
		CE SPACES 	//Composite element  	   	Acc. number and identifier for the corresponding composite element.
		composite_element
	) | (
		CX SPACES 	//Complex			a list of complexes which contain this factor
		complex
	) | (
		DE SPACES 	//Description 	  	short gene term (explicit gene name); GENE accession no.
		description
	) | (
		DR SPACES 	//External database links 	  	name of database (e. g. TRANSPRO, PathoDB, Flybase, EPD): database accession number; identifier (where available).
			//EMBL: accession number; identifier (first e position:last site position).
			//RSNP: accession number; EMBL: accession number; pos: SNP position in EMBL sequence; var: variation introduced by SNP.
		external_databases
	) | (
		EV SPACES //Evidence	Accession number of the evidence for this CE. Evidence for a given CE is an experiment of a certain type carried out with two individual interactions within a particular cell type.
		evidence
	) | (
		FA SPACES //Factor name   (normally the most commonly used) name of the transcription factor (NOTE: Greek letters are expanded to alpha, beta, gamma etc.) 
		factor_name
	) | (
		GE SPACES //GENE accession no.; short gene term; HGNC: standard gene symbol.
		gene
	) | (
		HC SPACES //Subfamilies  	   	lists entries, e.g. splice variants or family members of this isogroup/family entry
		sub_family
	) | (
		HP SPACES //Superfamilies  	   	lists generic entries (isogroup or family) to which this factor belongs
		super_family
	) | (
		ID SPACES //Identifier
		identifier
	) | (
		NA SPACES //Factor name  	   	designation of the binding transcription factor
		name
	) | (
		NUMBER_CODE SPACES 	//Binding Matrix  	 nucleotide frequency matrix with matrix head (A C G T) and derived IUPAC consensus in the last column
		number_code
	) | (
		MX SPACES //Matrices  	   	MATRIX table entries providing DNA-binding profiles of the factor (linked accession number; identifier.)
		matrices
	) | (
		OC SPACES   //Taxonomic classification  	   	systematic biological classification of the species
		taxonomy
	) | (
		PS SPACES 	//Position  	Identifies the position of the composite element within the promoter.
		position
	) | (
		PW SPACES 	//Pathways  	   	Indicates all the pathways and chains in which the respective molecule is involved.
		pathway
	) | (
		S1 SPACES 	//Reference point
		reference_point
	) | (
		SD SPACES 	//Start point
		short_description
	) | (
		SF SPACES 	//Start point
		start_position
	) | (
		SQ SPACES 	//Sequence 	  	Site sequence(s)
		sequence
	) | (
		ST SPACES 	//End point
		end_position
	) | (
		SY SPACES 	//alternative names of the transcription factor
		synonyms
	) | (
		TY SPACES 	//Type
		type
	)
	;

	
/** Override the following rules to parse the given fields. */
protected
accession_number : table_link_value[e->accession_number] ;

value_value [ std::string & s_arg ]
	: s:WORD { s_arg = s->getText(); } 
	  (
		(t:WS {s_arg += t -> getText();}
		 u:WORD {s_arg += u -> getText();})
 	   )*
 	   (WS)*
	;

binding_factor
{ FactorLinkPtr link(new FactorLink); std::string value; }
	:
	{ sps.push_lexer("list"); }
	{ sps.push_lexer("link"); }
	table_link_value [link->link]
	{ sps.pop_lexer(); }
	SEMI
	str_value[link->name]
	(
		SEMI
		{ sps.push_lexer("name_value_pair"); }
		(
			(WS)?
			(
				(QUALITY WS) { value_value(value); link->quality = atoi(value.c_str());}
				| (SPECIES WS) { value_value(value); link->species.push_back(value);}
				| (SITES WS INCLUDED WS) { value_value(value); link->sites_included = "yes" == value;}
				| (CELLULAR WS SOURCE WS) { value_value(value); link->cellular_source = value;}
			)
		)*
		{ sps.pop_lexer(); }
	)*
	{ e->factor_links.push_back(link); }
	{ sps.pop_lexer(); }
	;

protected
matrix_basis : ignored_value ;

protected
binding_site : ignored_value ;
	
protected
composite_element : ignored_value ;

protected
comment : ignored_value ;

protected
complex : ignored_value ;
	
protected
description
	:
	{ sps.push_lexer("string"); }
	s:STRING
	{ sps.pop_lexer(); }
	{ e->description = s->getText(); }
	;
	
protected
evidence : ignored_value ;
	
protected
factor_name : ignored_value ;
	
protected
gene : ignored_value ;
	
protected
name : ignored_value ;
	
protected
identifier : identifier_value[e->id] ;

protected
number_code : ignored_value ;
	
protected
matrices : ignored_value ;
	
protected
pathway : ignored_value ;

protected
position : ignored_value ;

protected
reference_point : ignored_value ;

protected
start_position : ignored_value ;

protected
end_position : ignored_value ;

protected
sequence : sequence_value[e->sequence] ;

protected
short_description : ignored_value ;

protected
super_family : ignored_value ;

protected
sub_family : ignored_value ;

protected
synonyms : ignored_value ;

protected
taxonomy : ignored_value ;

protected
type : ignored_value ;


/** Library of rules to use when overriding those above. */
protected
ignored_value
	: { sps.push_lexer("string"); } STRING { sps.pop_lexer(); } ; //parse a string and ignore it
	
identifier_value [ Identifier & id ]
	: { sps.push_lexer("id"); }
	str_value [id.species_group]
	DOLLAR (str_value [id.factor])?
	(
		DOLLAR
		str_value [id.factor]
	)?
	(UNDER { sps.pop_lexer(); sps.push_lexer("string"); } str_value [id.discriminating_extension])?
	{ sps.pop_lexer(); }
	;

table_link_value [ TableLink & tl ]
{ int n; }
	:
	{ sps.push_lexer("link"); }
	id:TABLE_ID { tl.table_id = parse_biobase_table(id->getText()); }
	{ sps.pop_lexer(); }
	n = integer_value { tl.entry_idx = n; }
	;

/** Appends the matched sequence to the given sequence returns whether a dot was matched (i.e. if the sequence was completed)
 */
sequence_value [ seq_t & seq ]
{
	bool bad_data = false;
	std::stringstream seq_stream;
	unsigned int start = 0; //initialise to avoid compiler warning
}
	:
	{ sps.push_lexer("sequence"); 
	  start = sps.markInputstream();}		//	So that we can go back if there is an error
	(
		(
			//match a sequence with spaces in it and an optional number or dot at the end.
			(
				s:SEQUENCE { seq_stream << s->getText(); }
				(
					SPACES
					|
					( //e.g. [1:30]
						GAP_START
						NUMBER
						COLON
						NUMBER
						GAP_END
						{ bad_data=true; } //we don't cater for sequences with gaps
					)
					|
					DOT
				)?
			)*
			(NUMBER)?
		)
		|
		(
			//empty sequence
			DOT { if (seq.size() == 0) bad_data=true; }
		)
	)
	{ 
		sps.pop_lexer(); 
		if ( bad_data)
			err_msg(start,"SQ");
		else
		{
			sps.commitInputstream();
			seq += seq_stream.str();
		}
	}
	;
	exception
	catch [const ANTLR_USE_NAMESPACE(antlr)ANTLRException & ex] 
	{   sps.pop_lexer(); 
		err_msg(start,"SQ");}
	
integer_value returns [ int i = -1 ]
	: { sps.push_lexer("number"); }
	n:NUMBER { i = sps.parse_number(n->getText()); }
	{ sps.pop_lexer(); } ;
	
str_value [ std::string & s_arg ]
	: s:STRING { s_arg = s->getText(); } ;

err_msg [ unsigned int mark, const char * field] 
	: { sps.rewindInputstream(mark);	//	Go back to where we started reading in the sequence
		inputState -> reset();			//  Remove any previously read tokens from the queue
		sps.push_lexer("string");}
		s:STRING 
	  {  std::cerr << "Problem parsing " << field << " field in " << e->get_name() << "\n" << 
				  field << " field = " << s -> getText() << "\n"; 
			 	 sps.pop_lexer();} ;
	
protected
external_databases
	:
	{
		sps.push_lexer( "colondotsemi" );
		static DatabaseRefParser parser( getInputState() );
		parser.database_ref();
		if( parser.parsed )
		{
			try
			{
				e->database_refs.push_back( parser.ref() );
			}
			catch ( std::exception const & e )
			{
				std::cerr << "Tried to parse \"" << parser.acc << "\" as " << parser.db << ": " << e.what() << "\n";
			}
		}
		sps.pop_lexer();
	}
	;
	
protected
super_family_value [ TableLinkVec & vec ]
{ TableLink l; }
	:
	{ sps.push_lexer("pathway"); }
	LESS
	{ sps.push_lexer("link"); }
	table_link_value[l] { vec.push_back(l); }
	{ sps.pop_lexer(); }
	GREATER
	{ sps.pop_lexer(); }
	{ sps.push_lexer("string"); }
	STRING
	{ sps.pop_lexer(); }
	;

