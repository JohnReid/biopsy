/* Copyright John Reid 2007
*/

#include "bio-pch.h"


#include "bio/defs.h"
#include "bio/database_ref.h"
#include "bio/common.h"

#include <boost/functional/hash.hpp>
#include <boost/algorithm/string/predicate.hpp>

BIO_NS_START

template< typename T, typename Str >
T
my_lexical_cast( const Str & s )
{
	try
	{
		return boost::lexical_cast< T >( s );
	}
	catch( const std::exception & e )
	{
		throw std::logic_error( BIO_MAKE_STRING( "\"" << s << "\": " << e.what() ) );
	}
}

const char * get_database_name( Database database )
{
	switch(database) {
	case AFFY_PROBE_DB: return "AFFY";
	case BKL_DB: return "BKL";
	case DIP_DB: return "DIP";
	case EMBL_DB: return "EMBL";
	case ENSEMBL_DB: return "ENSEMBL";
	case ENTREZ_GENE_DB: return "ENTREZ_GENE";
	case ENTREZ_PROTEIN_DB: return "ENTREZ_PROTEIN";
	case EPD_DB: return "EPD";
	case FLYBASE_DB: return "FLYBASE";
	case INPARANOID_DB: return "INPARANOID";
	case JASPAR_DB: return "JASPAR";
	case MGI_DB: return "MGI";
	case PATHO_DB: return "PATHO";
	case PDB_DB: return "PDB";
	case PIR_DB: return "PIR";
	case PROSITE_DB: return "PROSITE";
	case REFSEQ_DB: return "REFSEQ";
	case RGD_DB: return "RGD";
	case RSNP_DB: return "RSNP";
	case SGD_DB: return "SGD";
	case SMART_DB: return "SMART";
	case SWISSPROT_DB: return "SWISSPROT";
	case TAIR_DB: return "TAIR";
	case TRANSFAC_DB: return "TRANSFAC";
	case TRANSCOMPEL_DB: return "TRANSCOMPEL";
	case TRANSPATH_DB: return "TRANSPATH";
	case TRANSPRO_DB: return "TRANSPRO";
	case UNIGENE_DB: return "UNIGENE";
	case WORMBASE_DB: return "WORMBASE";
	case ZFIN_DB: return "ZFIN";
	case UNKNOWN_DB: return "<unknown db>";
	default: throw std::logic_error( BIO_MAKE_STRING( "Bad database enum: " << int( database ) ) );
	}
}

std::ostream &
operator<<( std::ostream & os, Database database )
{
	return os << get_database_name( database );
}


Database
parse_database_string( const std::string & database )
{
	using namespace boost::algorithm;
	if( iequals( database, "AFFYMETRIX" ) ) return AFFY_PROBE_DB;
	if( iequals( database, "BKL" ) ) return BKL_DB;
	if( iequals( database, "DIP" ) ) return DIP_DB;
	if( iequals( database, "EMBL" ) ) return EMBL_DB;
	if( iequals( database, "ENSEMBL" ) ) return ENSEMBL_DB;
	if( iequals( database, "ENTREZGENE" ) ) return ENTREZ_GENE_DB;
	if( iequals( database, "ENTREZPROTEIN" ) ) return ENTREZ_PROTEIN_DB;
	if( iequals( database, "EPD" ) ) return EPD_DB;
	if( iequals( database, "FLYBASE" ) || iequals( database, "FB" ) ) return FLYBASE_DB;
	if( iequals( database, "MGI" ) ) return MGI_DB;
	if( iequals( database, "PATHODB" ) ) return PATHO_DB;
	if( iequals( database, "PDB" ) ) return PDB_DB;
	if( iequals( database, "PIR" ) ) return PIR_DB;
	if( iequals( database, "REFSEQ" ) ) return REFSEQ_DB;
	if( iequals( database, "RGD" ) ) return RGD_DB;
	if( iequals( database, "RSNP" ) ) return RSNP_DB;
	if( iequals( database, "SMARTDB" ) ) return SMART_DB;
	if( iequals( database, "SWISSPROT" ) || iequals( database, "UNIPROT" ) ) return SWISSPROT_DB;
	if( iequals( database, "TRANSCOMPEL" ) ) return TRANSCOMPEL_DB;
	if( iequals( database, "TRANSFAC" ) ) return TRANSFAC_DB;
	if( iequals( database, "TRANSPATH" ) ) return TRANSPATH_DB;
	if( iequals( database, "UNIGENE" ) ) return UNIGENE_DB;
	if( iequals( database, "WB" ) ) return WORMBASE_DB;
	if( iequals( database, "ZFIN" ) ) return ZFIN_DB;
	return UNKNOWN_DB;
}


bool
db_ref::operator==(const db_ref & rhs) const
{
	return db == rhs.db && table == rhs.table && acc == rhs.acc;
}

bool
db_ref::operator<(const db_ref & rhs) const
{
	if( db < rhs.db ) return true;
	if( db > rhs.db ) return false;
	if( table < rhs.table ) return true;
	if( table > rhs.table ) return false;
	if( acc < rhs.acc ) return true;
	return false;
}


db_ref::db_ref()
: db( UNKNOWN_DB )
, acc( -1 )
{
}

db_ref::db_ref( 
	Database db,
	const std::string & table,
	int acc )
	: db( db )
	, table( table )
	, acc( acc )
{
}

int 
db_ref::compare( const db_ref & rhs ) const
{
	if( this->db <  rhs.db ) return -1;
	if( this->db != rhs.db ) return  1;
	if( this->table <  rhs.table ) return -1;
	if( this->table != rhs.table ) return  1;
	if( this->acc <  rhs.acc ) return -1;
	if( this->acc != rhs.acc ) return  1;
	return 0;
}


std::size_t hash_value( const db_ref & ref )
{
	std::size_t seed = 0;
	boost::hash_combine( seed, int( ref.db ) );
	boost::hash_combine( seed, ref.table );
	boost::hash_combine( seed, ref.acc );

	return seed;
}

/*
GenMAPP accepted ID's
Gene ID System 	System Code 	Applicable Species 	Example
Ensembl 	En 	B.taurus, C.elegans, C.familiaris, D.melanogaster, D.rerio, G.gallus, H.sapiens, M.musculus, R.norvegicus, S.cerevisiae 	ENSMUSG00000027793
UniProt/TrEMBL 	S 	B.taurus, C.elegans, C.familiaris, D.melanogaster, D.rerio, G.gallus, H.sapiens, M.musculus, R.norvegicus,  S.cerevisiae A2A2_HUMAN or O94973
Entrez Gene 	L 	B.taurus, C.elegans, C.familiaris, D.melanogaster, D.rerio, G.gallus, H.sapiens, M.musculus,  R.norvegicus,  S.cerevisiae 	68377
RefSeq (NM_xxxxxx only) 	Q 	B.taurus, C.elegans, C.familiaris, D.melanogaster, D.rerio, G.gallus, H.sapiens, M.musculus, R.norvegicus, S.cerevisiae 	NP_031407, NM_016749
Unigene 	U 	B.taurus, C.elegans, C.familiaris, D.melanogaster, D.rerio, G.gallus, H.sapiens, M.musculus, R.norvegicus 	Hs.451376
Affymetrix Probe Set ID (Affy) 	X 	B.taurus, C.elegans, C.familiaris, D.melanogaster, D.rerio, G.gallus, H.sapiens, M.musculus, R.norvegicus, S.cerevisiae 	100014_at
WormBase 	W 	C.elegans 	CE00005
FlyBase 	F 	D.melanogaster 	FBgn0000043
ZFIN 	Z 	D.rerio 	ZDB-GENE-000329-1
Mouse Genome Informatics (MGI) 	M 	M.musculus 	MGI:1194500
Rat Genome Database (RGD) 	R 	R.norvegicus 	RGD:70907
Saccharomyces Genome Database (SGD) 	D 	S.cerevisiae 	S0000157
PDB 	Pd 	B.taurus, C.elegans, C.familiaris, D.melanogaster, D.rerio, G.gallus, H.sapiens, M.musculus, R.norvegicus,  S.cerevisiae 	15C8
EMBL 	Em 	B.taurus, C.elegans, C.familiaris, D.melanogaster, D.rerio, G.gallus, H.sapiens, M.musculus,  R.norvegicus, S.cerevisiae 	AA000715
Other 	O 	- 	12345
*/

namespace detail {
	using namespace boost;
	static const regex affy_probe_acc_regex( " *([0-9]{6})_at *" );
	static const regex dip_acc_regex( " *DIP:([0-9A-Z]+) *" );
	static const regex embl_acc_regex( " *([A-Za-z_0-9\\.]+) *" );
	static const regex embl_acc_regex_2( // http://www.ncbi.nlm.nih.gov/Sequin/acc.html
		" *(?:"
		"(?:([A-Za-z]{1})([0-9]{5}))" //1 letter + 5 numerals
		"|(?:([A-Za-z]{2})([0-9]{6}))" //2 letters + 6 numerals
		"|(?:([A-Za-z]{3})([0-9]{5}))" //3 letters + 5 numerals
		"|(?:([A-Za-z]{4})([0-9]{8,10}))" //4 letters + 8-10 numerals
		") *"
		);
	static const regex ensembl_acc_regex( " *([A-Z]{2,})([0-9]{3,14}) *" );
	static const regex entrez_gene_acc_regex( " *([0-9]+) *" );
	static const regex entrez_protein_acc_regex( " *([0-9]+) *" );
	static const regex flybase_acc_regex( " *FBgn([0-9]{7}) *" );
	static const regex mgi_acc_regex( " *(?:MGI:)?((?:[0-9_\\(\\)\\.]|-)+) *" );
	static const regex pdb_acc_regex( " *([0-9]{2}[A-Z][0-9]) *" );
	static const regex refseq_acc_regex( " *(N[PM])_([0-9]{6}) *" ); //http://www.ncbi.nlm.nih.gov/RefSeq/key.html#accession
	static const regex rgd_acc_regex( " *RGD:([0-9]+) *" );
	static const regex sgd_acc_regex( " *S([0-9]{7}) *" );
	static const regex swissprot_acc_regex( " *([A-Z0-9]+)(?:-([0-9]+))? *" );
	static const regex swissprot_acc_regex_2( " *((?:(?:[A-NR-Z][0-9][A-Z])|(?:[OPQ][0-9][A-Z0-9]))[A-Z0-9][A-Z0-9][0-9])(?:-([0-9]+))? *" ); //see http://expasy.org/sprot/userman.html#AC_line
	static const regex trans_acc_regex( " *([A-Z]+)([0-9]+) *" );
	static const regex unigene_acc_regex( " *([A-Z][a-z]).([0-9]{6}) *" );
	static const regex wormbase_acc_regex( " *CE([0-9]{6}) *" );
	static const regex zfin_acc_regex( " *ZDB-GENE-([0-9]{6})-([0-9]+) *" );
} // namespace detail

db_ref 
parse_db_ref_as( const std::string & acc, Database db )
{
	using namespace boost;
	cmatch what;
	switch( db )
	{
	case DIP_DB:
		if( ! regex_match( acc.c_str(), what, detail::dip_acc_regex ) ) break;
		else return db_ref( db, what[1], -1 );

	case EMBL_DB:
		if( ! regex_match( acc.c_str(), what, detail::embl_acc_regex_2 ) ) break;
		else return db_ref( db, what[1], my_lexical_cast< int >( what[2] ) );

	case ENSEMBL_DB:
		if( ! regex_match( acc.c_str(), what, detail::ensembl_acc_regex ) ) break;
		else return db_ref( db, what[1], my_lexical_cast< int >( what[2] ) );

	case ENTREZ_GENE_DB:
		if( ! regex_match( acc.c_str(), what, detail::entrez_gene_acc_regex ) ) break;
		else return db_ref( db, "", my_lexical_cast< int >( what[1] ) );

	case ENTREZ_PROTEIN_DB:
		if( ! regex_match( acc.c_str(), what, detail::entrez_protein_acc_regex ) ) break;
		else return db_ref( db, "", my_lexical_cast< int >( what[1] ) );

	case FLYBASE_DB:
		if( ! regex_match( acc.c_str(), what, detail::flybase_acc_regex ) ) break;
		else return db_ref( db, "", my_lexical_cast< int >( what[1] ) );

	case MGI_DB:
		if( ! regex_match( acc.c_str(), what, detail::mgi_acc_regex ) ) break;
		else return db_ref( db, "", my_lexical_cast< int >( what[1] ) );

	case SWISSPROT_DB:
		if( ! regex_match( acc.c_str(), what, detail::swissprot_acc_regex_2 ) ) break;
		else return db_ref( db, what[1], -1 );

	case TRANSCOMPEL_DB:
	case TRANSFAC_DB:
	case TRANSPATH_DB:
		if( ! regex_match( acc.c_str(), what, detail::trans_acc_regex ) ) break;
		else return db_ref( db, what[1], my_lexical_cast< int >( what[2] ) );

	default:
		throw std::logic_error( BIO_MAKE_STRING( "Don't know how to parse database: " << db ) );
	}
	throw std::logic_error( BIO_MAKE_STRING( "Could not parse accession \"" << acc << "\" for database: " << db ) );
}



std::ostream &
operator<<( std::ostream & os, const db_ref & ref )
{
	boost::io::ios_fill_saver ifs( os );

	switch( ref.db )
	{
	case AFFY_PROBE_DB: return os << std::setw( 6 ) << std::setfill( '0' ) << ref.acc << "_at";
	case DIP_DB: return os << "DIP:" << ref.table;
	case EMBL_DB: return os << "EMBL:" << ref.table << std::setw( 6 ) << std::setfill( '0' ) << ref.acc;
	case ENSEMBL_DB: return os << ref.table << std::setw(11) << std::setfill( '0' ) << ref.acc;
	case ENTREZ_GENE_DB: return os << "ENTREZ_GENE:" << std::setw( 5 ) << std::setfill( '0' ) << ref.acc;
	case ENTREZ_PROTEIN_DB: return os << "ENTREZ_PROTEIN:" << std::setw( 5 ) << std::setfill( '0' ) << ref.acc;
	case FLYBASE_DB: return os << "FBgn" << std::setw( 7 ) << std::setfill( '0' ) << ref.acc;
	case MGI_DB: return os << "MGI:" << ref.acc;
	case PDB_DB: return os << "PDB:" << ref.table;
	case REFSEQ_DB: return os << "REFSEQ:" << ref.table << std::setw( 6 ) << std::setfill( '0' ) << ref.acc;
	case RGD_DB: return os << "RGD:" << ref.acc;
	case SGD_DB: return os << "SGD:" << ref.table;
	case SWISSPROT_DB: return os << "SWISSPROT:" << ref.table;
	case UNIGENE_DB: return os << "UNIGENE:" << ref.table << "." << std::setw( 6 ) << std::setfill( '0' ) << ref.acc;
	case WORMBASE_DB: return os << "CE" << std::setw( 6 ) << std::setfill( '0' ) << ref.acc;
	case ZFIN_DB: return os << "ZDB-GENE-" << ref.acc << "-" << ref.table;
	default:
		return os << ref.db << ":" << ref.table << ":" << ref.acc;
	}
}

db_ref
try_to_parse_db_ref( const std::string & acc )
{
	using namespace boost;

	try
	{
		cmatch what; 

		//ensembl
		if( regex_match( acc.c_str(), what, detail::ensembl_acc_regex ) )
		{
			return 
				db_ref( 
					ENSEMBL_DB, 
					what[1],
					my_lexical_cast< int >( what[2] ) );
		}

		//swissprot 1
		if( regex_match( acc.c_str(), what, detail::swissprot_acc_regex_2 ) )
		{
			return 
				db_ref( 
					SWISSPROT_DB, 
					what[1],
					0 );
		}

		//entrez_gene
		if( regex_match( acc.c_str(), what, detail::entrez_gene_acc_regex ) )
		{
			return 
				db_ref( 
					ENTREZ_GENE_DB, 
					"",
					my_lexical_cast< int >( what[1] ) );
		}

		//entrez_protein
		if( regex_match( acc.c_str(), what, detail::entrez_protein_acc_regex ) )
		{
			return 
				db_ref( 
					ENTREZ_PROTEIN_DB, 
					"",
					my_lexical_cast< int >( what[1] ) );
		}

		//refseq
		if( regex_match( acc.c_str(), what, detail::refseq_acc_regex ) )
		{
			return 
				db_ref( 
					REFSEQ_DB, 
					what[1],
					my_lexical_cast< int >( what[2] ) );
		}

		//unigene
		if( regex_match( acc.c_str(), what, detail::unigene_acc_regex ) )
		{
			return 
				db_ref( 
					UNIGENE_DB, 
					what[1],
					my_lexical_cast< int >( what[2] ) );
		}

		//affy_probe
		if( regex_match( acc.c_str(), what, detail::affy_probe_acc_regex ) )
		{
			return 
				db_ref( 
					AFFY_PROBE_DB, 
					"",
					my_lexical_cast< int >( what[1] ) );
		}

		//wormbase
		if( regex_match( acc.c_str(), what, detail::wormbase_acc_regex ) )
		{
			return 
				db_ref( 
					WORMBASE_DB, 
					"",
					my_lexical_cast< int >( what[1] ) );
		}

		//flybase
		if( regex_match( acc.c_str(), what, detail::flybase_acc_regex ) )
		{
			return 
				db_ref( 
					FLYBASE_DB, 
					"",
					my_lexical_cast< int >( what[1] ) );
		}

		//zfin
		if( regex_match( acc.c_str(), what, detail::zfin_acc_regex ) )
		{
			return 
				db_ref( 
					ZFIN_DB, 
					what[2],
					my_lexical_cast< int >( what[1] ) );
		}

		//mgi
		if( regex_match( acc.c_str(), what, detail::mgi_acc_regex ) )
		{
			return 
				db_ref( 
					MGI_DB, 
					"",
					my_lexical_cast< int >( what[1] ) );
		}

		//rgd
		if( regex_match( acc.c_str(), what, detail::rgd_acc_regex ) )
		{
			return 
				db_ref( 
					RGD_DB, 
					"",
					my_lexical_cast< int >( what[1] ) );
		}

		//pdb
		if( regex_match( acc.c_str(), what, detail::pdb_acc_regex ) )
		{
			return 
				db_ref( 
					PDB_DB, 
					what[1],
					0 );
		}

		//embl
		if( regex_match( acc.c_str(), what, detail::embl_acc_regex_2 ) )
		{
			return 
				db_ref( 
					EMBL_DB, 
					what[1],
					my_lexical_cast< int >( what[2] ) );
		}

		//dip
		if( regex_match( acc.c_str(), what, detail::dip_acc_regex ) )
		{
			return 
				db_ref( 
					DIP_DB, 
					what[1],
					-1 );
		}

	}
	catch( ... )
	{
	}

	return db_ref( UNKNOWN_DB, acc, 0 );
}

db_ref 
parse_db_ref( const std::string & acc )
{
	db_ref result = try_to_parse_db_ref( acc );

	if( UNKNOWN_DB == result.db )
	{
		throw std::logic_error( BIO_MAKE_STRING( "Could not parse accession: \"" << acc << "\"" ) );
	}

	return result;
}


db_ref
db_ref_from_transfac_table_link( const BIO_NS::TableLink & link )
{
	return db_ref(
		database_from_trans_data_type( link.table_id ),
		BIO_MAKE_STRING( link.table_id ),
		link.entry_idx );
}

TableLink
transfac_table_link_from_db_ref( const db_ref & ref )
{
	if( TRANSFAC_DB == ref.db
		|| TRANSCOMPEL_DB == ref.db
		|| TRANSPATH_DB == ref.db )
	{
		return TableLink( parse_biobase_table( ref.table ), ref.acc );
	}
	else
	{
		return TableLink();
	}
}

const char * database_name_for_ensembl_table( const std::string & table )
{
	using namespace boost::algorithm;
	if( iequals( table, "ENSMUSG" ) ) return "Mus_musculus";
	if( iequals( table, "ENSDARG" ) ) return "Danio_rerio";
	if( iequals( table, "CG" ) ) return "Drosophila_melanogaster";
	if( iequals( table, "ENSGALG" ) ) return "Gallus_gallus";
	if( iequals( table, "ENSG" ) ) return "Homo_sapiens";
	if( iequals( table, "ENSMMUG" ) ) return "Macaca_mulatta";
	if( iequals( table, "ENSRNOG" ) ) return "Rattus_norvegicus";
	if( iequals( table, "SINFRUG" ) ) return "Takifugu_rubripes";
	if( iequals( table, "GSTENG" ) ) return "Tetraodon_nigroviridis";
	if( iequals( table, "ENSXETG" ) ) return "Xenopus_tropicalis";
	throw std::logic_error( BIO_MAKE_STRING( "Do not know database name for Ensembl table: " << table ) );
}

std::string url_for( const db_ref & ref )
{
	switch( ref.db )
	{
	case SWISSPROT_DB: //http://expasy.org/uniprot/P22557
		return BIO_MAKE_STRING( "http://expasy.org/uniprot/" << ref.table ) ;

	case ENTREZ_PROTEIN_DB: //http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?db=protein&val=33466082
		return BIO_MAKE_STRING( "http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?db=protein&val=" << ref.acc );

	case ENTREZ_GENE_DB: //http://www.ncbi.nlm.nih.gov/sites/entrez?Db=gene&Cmd=ShowDetailView&TermToSearch=72568
		return BIO_MAKE_STRING( "http://www.ncbi.nlm.nih.gov/sites/entrez?Db=gene&Cmd=ShowDetailView&TermToSearch=" << ref.acc );

	case ENSEMBL_DB: //http://www.ensembl.org/Mus_musculus/geneview?gene=ENSMUSG00000008976
		return BIO_MAKE_STRING( 
			"http://www.ensembl.org/" 
			<< database_name_for_ensembl_table(ref.table)
			<< "/geneview?gene=" << ref );

	case FLYBASE_DB: //http://www.flybase.net/reports/FBgn0000043.html
		return BIO_MAKE_STRING( "http://www.flybase.net/reports/" << ref << ".html" );

	case MGI_DB: //http://www.informatics.jax.org/searches/accession_report.cgi?id=MGI:97275
		return BIO_MAKE_STRING( "http://www.informatics.jax.org/searches/accession_report.cgi?id=" << ref );

	case TRANSCOMPEL_DB:
	case TRANSFAC_DB:
	case TRANSPATH_DB: //http://www.biobase-international.com/cgi-bin/biobase/transpath/7.4/bin/start.cgi?ac=MO000026428
		return transfac_table_link_from_db_ref( ref ).get_url();

	default:
		throw std::logic_error( BIO_MAKE_STRING( "Do not know url for this database: " << ref.db ) );
	}
}


std::pair< bool, db_ref >
parse_transfac_db_ref( const std::string & ref_string )
{
	using namespace boost;
	const static regex ref_regex( "([A-Za-z0-9/-]+):([A-Za-z0-9_\\-:/\\.]+)" );

	cmatch what;
	if( ! regex_match( ref_string.c_str(), what, ref_regex ) )
	{
		//throw std::logic_error( BIO_MAKE_STRING( "Could not parse ref: \"" << ref_string << "\"" ) );
		return std::make_pair( false, db_ref() );
	}
	else
	{
		const Database db = parse_database_string( what[1] );
		switch( db )
		{
		case UNKNOWN_DB:
			return std::make_pair( false, db_ref( db, "", -1 ) );
		case AFFY_PROBE_DB:
		case BKL_DB:
		case PDB_DB:
		case PIR_DB:
		case REFSEQ_DB:
		case RGD_DB:
		case UNIGENE_DB:
			return std::make_pair( false, db_ref( db, "", -1 ) );

		default:
			return std::make_pair( true, parse_db_ref_as( what[2], db ) );
		}
	}
}

std::pair< bool, db_ref >
parse_transfac_db_ref_line( const std::string & line )
{
	using namespace boost;
	const static std::string ref_regex_string = "[A-Za-z0-9_\\-:/\\.]+";
	const static regex molecule_dr_regex_1( BIO_MAKE_STRING( "<(" << ref_regex_string << ")(?:>)?" ) );
	const static regex molecule_dr_regex_2( BIO_MAKE_STRING( "\\{[A-Z]+\\}(" << ref_regex_string << ")" ) );

	cmatch what;
	std::string ref_string;
	if( regex_search( line.c_str(), what, molecule_dr_regex_1 ) )
	{
		ref_string = what[1];
	}
	else if ( regex_search( line.c_str(), what, molecule_dr_regex_2 ) )
	{
		ref_string = what[1];
	}
	else
	{
		//cout << "Could not parse: \"" << line << "\"\n";
		return std::make_pair( false, db_ref() );
	}

	return parse_transfac_db_ref( ref_string );
}



BIO_NS_END

