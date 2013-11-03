/**
@file

Copyright John Reid 2006

*/

#include "bio-pch.h"

//#define BOOST_SPIRIT_DEBUG

#include <bio/spirit/transfac_qi.h>

#include <boost/assign.hpp>

#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>

/** Test a parser. */
template< typename P >
void
test_parser(
    char const * input,
    P const & p,
    bool full_match = true
) {
    char const * f( input );
    char const * l( f + strlen( f ) );
    BOOST_CHECK_MESSAGE( ::boost::spirit::qi::parse( f, l, p ), "Problem parsing \"" << input << "\""  );
    BOOST_CHECK_MESSAGE(
    	! full_match || f == l,
    	"Did not fully match \"" << input << "\""
    	<< " only matched \"" << std::string( input, f ) << "\""
    );
}

/** Test a phrase parser. */
template< typename P >
void
test_phrase_parser(
    char const * input,
    P const & p,
    bool full_match = true
) {
	using namespace ::boost::spirit::qi;
    char const * f( input );
    char const * l( f + strlen( f ) );
    BOOST_CHECK_MESSAGE( phrase_parse( f, l, p, ascii::space ), "Problem parsing \"" << input << "\"" );
    BOOST_CHECK_MESSAGE(
    	! full_match || f == l,
    	"Did not fully match \"" << input << "\""
    	<< " only matched \"" << std::string( input, f ) << "\""
    );
}



/** Test a parser. */
template< typename P, typename T >
T &
test_parser_attr(
    char const * input,
    P const & p,
    T & attr,
    bool full_match = true
) {
    char const * f( input );
    char const * l( f + strlen( f ) );
    BOOST_CHECK_MESSAGE( ::boost::spirit::qi::parse( f, l, p, attr ), "Problem parsing \"" << input << "\""  );
    BOOST_CHECK_MESSAGE(
    	! full_match || f == l,
    	"Did not fully match \"" << input << "\""
    	<< " only matched \"" << std::string( input, f ) << "\""
    );
    return attr;
}

/** Test a phrase parser. */
template< typename P, typename T >
T &
test_phrase_parser_attr(
    char const * input,
    P const & p,
    T & attr,
    bool full_match = true
) {
	using namespace ::boost::spirit::qi;
    char const * f( input );
    char const * l( f + strlen( f ) );
    BOOST_CHECK_MESSAGE( phrase_parse( f, l, p, ascii::blank, attr ), "Problem parsing \"" << input << "\"" );
    BOOST_CHECK_MESSAGE(
    	! full_match || f == l,
    	"Did not fully match \"" << input << "\""
    	<< " only matched \"" << std::string( input, f ) << "\""
    );
    return attr;
}



BOOST_AUTO_TEST_CASE( ParseTransData )
{
	using namespace BIO_NS;
	BIO_NS::spirit::transdata_parser< char const * > parser;
	BIO_NS::TransData attr;

	BOOST_CHECK_EQUAL( test_phrase_parser_attr( "Cell", parser, attr ), CELL_DATA );
	BOOST_CHECK_EQUAL( test_phrase_parser_attr( "C", parser, attr ), COMPEL_DATA );
	BOOST_CHECK_EQUAL( test_phrase_parser_attr( "ev", parser, attr ), EVIDENCE_DATA );
	BOOST_CHECK_EQUAL( test_phrase_parser_attr( "T", parser, attr ), FACTOR_DATA );
	BOOST_CHECK_EQUAL( test_phrase_parser_attr( "FR", parser, attr ), FRAGMENT_DATA );
	BOOST_CHECK_EQUAL( test_phrase_parser_attr( "G", parser, attr ), GENE_DATA );
	BOOST_CHECK_EQUAL( test_phrase_parser_attr( "M", parser, attr ), MATRIX_DATA );
	BOOST_CHECK_EQUAL( test_phrase_parser_attr( "MO", parser, attr ), MOLECULE_DATA );
	BOOST_CHECK_EQUAL( test_phrase_parser_attr( "CH", parser, attr ), PATHWAY_DATA );
	BOOST_CHECK_EQUAL( test_phrase_parser_attr( "XN", parser, attr ), REACTION_DATA );
	BOOST_CHECK_EQUAL( test_phrase_parser_attr( "RE", parser, attr ), REFERENCE_DATA );
	BOOST_CHECK_EQUAL( test_phrase_parser_attr( "S", parser, attr ), S_DATA );
	BOOST_CHECK_EQUAL( test_phrase_parser_attr( "R", parser, attr ), SITE_DATA );
}

BOOST_AUTO_TEST_CASE( ParseTableLink )
{
	using namespace BIO_NS;
	BIO_NS::spirit::tablelink_parser< char const * > parser;
	BIO_NS::TableLink attr;

	BOOST_CHECK_EQUAL( test_phrase_parser_attr( "Cell01", parser, attr ), TableLink( CELL_DATA, 1 ) );
	BOOST_CHECK_EQUAL( test_phrase_parser_attr( "Compel1", parser, attr ), TableLink( COMPEL_DATA, 1 ) );
	BOOST_CHECK_EQUAL( test_phrase_parser_attr( "Evidence001", parser, attr ), TableLink( EVIDENCE_DATA, 1 ) );
	BOOST_CHECK_EQUAL( test_phrase_parser_attr( "Factor001", parser, attr ), TableLink( FACTOR_DATA, 1 ) );
	BOOST_CHECK_EQUAL( test_phrase_parser_attr( "Fragment001", parser, attr ), TableLink( FRAGMENT_DATA, 1 ) );
	BOOST_CHECK_EQUAL( test_phrase_parser_attr( "Gene001", parser, attr ), TableLink( GENE_DATA, 1 ) );
	BOOST_CHECK_EQUAL( test_phrase_parser_attr( "Matrix001", parser, attr ), TableLink( MATRIX_DATA, 1 ) );
	BOOST_CHECK_EQUAL( test_phrase_parser_attr( "Molecule001", parser, attr ), TableLink( MOLECULE_DATA, 1 ) );
	BOOST_CHECK_EQUAL( test_phrase_parser_attr( "Pathway001", parser, attr ), TableLink( PATHWAY_DATA, 1 ) );
	BOOST_CHECK_EQUAL( test_phrase_parser_attr( "Reaction001", parser, attr ), TableLink( REACTION_DATA, 1 ) );
	BOOST_CHECK_EQUAL( test_phrase_parser_attr( "Reference001", parser, attr ), TableLink( REFERENCE_DATA, 1 ) );
	BOOST_CHECK_EQUAL( test_phrase_parser_attr( "s001", parser, attr ), TableLink( S_DATA, 1 ) );
	BOOST_CHECK_EQUAL( test_phrase_parser_attr( "Site1", parser, attr ), TableLink( SITE_DATA, 1 ) );
}

BOOST_AUTO_TEST_CASE( ParseDatabase )
{
	using namespace BIO_NS;
	BIO_NS::spirit::db_parser< char const * > parser;
	BIO_NS::Database attr;

	BOOST_CHECK_EQUAL( test_phrase_parser_attr( "AFFYMETRIX", parser, attr ), AFFY_PROBE_DB );
	BOOST_CHECK_EQUAL( test_phrase_parser_attr( "BKL", parser, attr ), BKL_DB );
	BOOST_CHECK_EQUAL( test_phrase_parser_attr( "DIP", parser, attr ), DIP_DB );
	BOOST_CHECK_EQUAL( test_phrase_parser_attr( "EMBL", parser, attr ), EMBL_DB );
	BOOST_CHECK_EQUAL( test_phrase_parser_attr( "ENSEMBL", parser, attr ), ENSEMBL_DB );
	BOOST_CHECK_EQUAL( test_phrase_parser_attr( "ENTREZGENE", parser, attr ), ENTREZ_GENE_DB );
	BOOST_CHECK_EQUAL( test_phrase_parser_attr( "ENTREZPROTEIN", parser, attr ), ENTREZ_PROTEIN_DB );
	BOOST_CHECK_EQUAL( test_phrase_parser_attr( "EPD", parser, attr ), EPD_DB );
	BOOST_CHECK_EQUAL( test_phrase_parser_attr( "FLYBASE", parser, attr ), FLYBASE_DB );
	BOOST_CHECK_EQUAL( test_phrase_parser_attr( "MGI", parser, attr ), MGI_DB );
	BOOST_CHECK_EQUAL( test_phrase_parser_attr( "PATHODB", parser, attr ), PATHO_DB );
	BOOST_CHECK_EQUAL( test_phrase_parser_attr( "PDB", parser, attr ), PDB_DB );
	BOOST_CHECK_EQUAL( test_phrase_parser_attr( "PIR", parser, attr ), PIR_DB );
	BOOST_CHECK_EQUAL( test_phrase_parser_attr( "REFSEQ", parser, attr ), REFSEQ_DB );
	BOOST_CHECK_EQUAL( test_phrase_parser_attr( "RGD", parser, attr ), RGD_DB );
	BOOST_CHECK_EQUAL( test_phrase_parser_attr( "RSNP", parser, attr ), RSNP_DB );
	BOOST_CHECK_EQUAL( test_phrase_parser_attr( "SMARTDB", parser, attr ), SMART_DB );
	BOOST_CHECK_EQUAL( test_phrase_parser_attr( "SWISSPROT", parser, attr ), SWISSPROT_DB );
	BOOST_CHECK_EQUAL( test_phrase_parser_attr( "TRANSCOMPEL", parser, attr ), TRANSCOMPEL_DB );
	BOOST_CHECK_EQUAL( test_phrase_parser_attr( "TRANSFAC", parser, attr ), TRANSFAC_DB );
	BOOST_CHECK_EQUAL( test_phrase_parser_attr( "TRANSPATH", parser, attr ), TRANSPATH_DB );
	BOOST_CHECK_EQUAL( test_phrase_parser_attr( "UNIGENE", parser, attr ), UNIGENE_DB );
	BOOST_CHECK_EQUAL( test_phrase_parser_attr( "WB", parser, attr ), WORMBASE_DB );
	BOOST_CHECK_EQUAL( test_phrase_parser_attr( "ZFIN", parser, attr ), ZFIN_DB );
}

BOOST_AUTO_TEST_CASE( ParseDbs )
{
	using namespace BIO_NS;
	using namespace BIO_NS::spirit;
	db_ref attr;

	BOOST_CHECK_EQUAL(
		test_phrase_parser_attr( "M37787", embl_parser< char const * >(), attr=db_ref() ),
		db_ref( EMBL_DB, "M37787", 0 )
	);
	BOOST_CHECK_EQUAL(
		test_phrase_parser_attr( "X02996", embl_parser< char const * >(), attr=db_ref() ),
		db_ref( EMBL_DB, "X02996", 0 )
	);
	BOOST_CHECK_EQUAL(
		test_phrase_parser_attr( "AL137852.15", embl_parser< char const * >(), attr=db_ref() ),
		db_ref( EMBL_DB, "AL137852.15", 0 )
	);
	BOOST_CHECK_EQUAL(
		test_phrase_parser_attr( "DIP:(AC987987234)", dip_parser< char const * >(), attr=db_ref() ),
		db_ref( DIP_DB, "AC987987234", 0 )
	);
	BOOST_CHECK_EQUAL(
		test_phrase_parser_attr( "Z746_15_at", affyprobe_parser< char const * >(), attr=db_ref() ),
		db_ref( AFFY_PROBE_DB, "Z746_15_at", 0 )
	);
	BOOST_CHECK_EQUAL(
		test_phrase_parser_attr( "200924_3p_s_at.", affyprobe_parser< char const * >(), attr=db_ref() ),
		db_ref( AFFY_PROBE_DB, "200924_3p_s_at", 0 )
	);
	BOOST_CHECK_EQUAL(
		test_phrase_parser_attr( "Hs.172928.0.A1_3p_a_at", affyprobe_parser< char const * >(), attr=db_ref() ),
		db_ref( AFFY_PROBE_DB, "Hs.172928.0.A1_3p_a_at", 0 )
	);
}

BOOST_AUTO_TEST_CASE( ParseDbRef )
{
	using namespace BIO_NS;
	spirit::dbref_parser< char const * > parser;
	db_ref attr;

	BOOST_CHECK_EQUAL( test_phrase_parser_attr( "EMBL: J01901; XX2.", parser, attr=db_ref() ), db_ref( EMBL_DB, "J01901", 0 ) );
	BOOST_CHECK_EQUAL( test_phrase_parser_attr( "EMBL: AC008894.9 (70005:70792).", parser, attr=db_ref() ), db_ref( EMBL_DB, "AC008894.9", 0 ) );
	BOOST_CHECK_EQUAL( test_phrase_parser_attr( "ENTREZGENE: 396450.", parser, attr=db_ref() ), db_ref( ENTREZ_GENE_DB, "", 396450 ) );
	BOOST_CHECK_EQUAL( test_phrase_parser_attr( "ENTREZPROTEIN: 2.", parser, attr=db_ref() ), db_ref( ENTREZ_PROTEIN_DB, "", 2 ) );
	BOOST_CHECK_EQUAL( test_phrase_parser_attr( "FLYBASE: FBgn0002561; l(1)sc.", parser, attr=db_ref() ), db_ref( FLYBASE_DB, "", 2561 ) );
	BOOST_CHECK_EQUAL( test_phrase_parser_attr( "JASPAR: MA0100.", parser, attr=db_ref() ), db_ref( JASPAR_DB, "MA", 100 ) );
	BOOST_CHECK_EQUAL( test_phrase_parser_attr( "PROSITE: PS00035, PS00465, PS00027.", parser, attr=db_ref() ), db_ref( PROSITE_DB, "", 35 ) );
	BOOST_CHECK_EQUAL( test_phrase_parser_attr( "PROSITE: PS00038.", parser, attr=db_ref() ), db_ref( PROSITE_DB, "", 38 ) );
	BOOST_CHECK_EQUAL( test_phrase_parser_attr( "REFSEQ: XM_002343739.", parser, attr=db_ref() ), db_ref( REFSEQ_DB, "XM", 2343739 ) );
	BOOST_CHECK_EQUAL( test_phrase_parser_attr( "SWISSPROT: P10084.", parser, attr=db_ref() ), db_ref( SWISSPROT_DB, "P10084", 0 ) );
	BOOST_CHECK_EQUAL( test_phrase_parser_attr( "TRANSCOMPEL: C00024.", parser, attr=db_ref() ), db_ref( TRANSCOMPEL_DB, "C", 24 ) );
	BOOST_CHECK_EQUAL( test_phrase_parser_attr( "TRANSPATH: G000005.", parser, attr=db_ref() ), db_ref( TRANSPATH_DB, "G", 5 ) );
	BOOST_CHECK_EQUAL( test_phrase_parser_attr( "TRANSPRO: HSA_1699.", parser, attr=db_ref() ), db_ref( TRANSPRO_DB, "HSA_1699", 0 ) );
}

BOOST_AUTO_TEST_CASE( ParseIdentifier )
{
	using namespace BIO_NS;
	BIO_NS::spirit::identifier_parser< char const * > parser;

	{
		BIO_NS::Identifier attr;
		test_phrase_parser_attr( "RAT$MCK_02", parser, attr );
		BOOST_CHECK_EQUAL( attr.species_group, "RAT" );
		BOOST_CHECK_EQUAL( attr.factor, "MCK" );
		BOOST_CHECK_EQUAL( attr.discriminating_extension, "02" );
	}

	{
		BIO_NS::Identifier attr;
		test_phrase_parser_attr( "AAV$P5", parser, attr );
		BOOST_CHECK_EQUAL( attr.species_group, "AAV" );
		BOOST_CHECK_EQUAL( attr.factor, "P5" );
		BOOST_CHECK_EQUAL( attr.discriminating_extension, "" );
	}

	{
		BIO_NS::Identifier attr;
		test_phrase_parser_attr( "AS$_NFKAPPAB_46", parser, attr );
		BOOST_CHECK_EQUAL( attr.species_group, "AS" );
		BOOST_CHECK_EQUAL( attr.factor, "NFKAPPAB" );
		BOOST_CHECK_EQUAL( attr.discriminating_extension, "46" );
	}

	{
		BIO_NS::Identifier attr;
		test_phrase_parser_attr( "Cp-APRE", parser, attr );
		BOOST_CHECK_EQUAL( attr.species_group, "" );
		BOOST_CHECK_EQUAL( attr.factor, "Cp-APRE" );
		BOOST_CHECK_EQUAL( attr.discriminating_extension, "" );
	}

	{
		BIO_NS::Identifier attr;
		test_phrase_parser_attr( "HS$IFI6_01", parser, attr );
		BOOST_CHECK_EQUAL( attr.species_group, "HS" );
		BOOST_CHECK_EQUAL( attr.factor, "IFI6" );
		BOOST_CHECK_EQUAL( attr.discriminating_extension, "01" );
	}

	{
		BIO_NS::Identifier attr;
		test_phrase_parser_attr( "AS$RUNX3_", parser, attr );
		BOOST_CHECK_EQUAL( attr.species_group, "AS" );
		BOOST_CHECK_EQUAL( attr.factor, "RUNX3" );
		BOOST_CHECK_EQUAL( attr.discriminating_extension, "" );
	}

}

BOOST_AUTO_TEST_CASE( ParseFactorLink )
{
	using namespace BIO_NS;
	BIO_NS::spirit::factorlink_parser< char const * > parser;
	using namespace boost::assign;

	std::vector< char const * > test_inputs;
	test_inputs +=
		"T02223; AML2; Quality: 2; Species: human, Homo sapiens.",
		"T00140; c-Myc-isoform1; Quality: 6; Species: human, Homo sapiens; Cellular source: 0123, Jurkat.",
		"T00428; ISGF-3; Quality: 6; Species: human, Homo sapiens.",
		"T01091; CPRF-1; Species: parsley, Petroselinum crispum (P. hortense); site(s) included: yes.",
		"T01228; Sp1; Species: hamster, Cricetulus sp.; site(s) included: no.",
		"T01128; MyoD; Species: chick, Gallus gallus; site(s) included: no.",
		"T15547; MyoD; Species: chick, Gallus gallus; site(s) included: no.",
		"T00524; MyoD; Species: clawed frog, Xenopus laevis; site(s) included: no.",
		"T10010; MyoD; Species: clawed frog, Xenopus laevis; site(s) included: no.",
		"T03902; MyoD; Species: common carp, Cyprinus carpio; site(s) included: no.",
		"T00525; MyoD; Species: human, Homo sapiens; site(s) included: no.",
		"T09197; MyoD; Species: human, Homo sapiens; site(s) included: no.",
		"T10027; MyoD; Species: Mammalia; site(s) included: no.",
		"T00527; MyoD; Species: monkey, Cercopithecus aethiops; site(s) included: no.",
		"T00526; MyoD; Species: mouse, Mus musculus; site(s) included: yes.",
		"T09177; MyoD; Species: mouse, Mus musculus; site(s) included: yes.",
		"T01551; MyoD; Species: quail, Coturnix coturnix; site(s) included: no.",
		"T15571; MyoD; Species: rat, Rattus norvegicus; site(s) included: no.",
		"T15572; MyoD; Species: rat, Rattus norvegicus; site(s) included: no.",
		"T03907; MyoD; Species: zebra fish, Brachydanio rerio; site(s) included: no.",
		"T03912; MyoD (275 AA); Species: rainbow trout, Oncorhynchus mykiss; site(s) included: no.",
		"T03911; MyoD (376 AA); Species: rainbow trout, Oncorhynchus mykiss; site(s) included: no."
		;

	BOOST_FOREACH( char const * input, test_inputs ) {
		BIO_NS::FactorLink attr;
		test_phrase_parser_attr( input, parser, attr );
		//BIO_NS::spirit::print_factorlink( attr );
	}
}
