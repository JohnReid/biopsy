/** Copyright John Reid 2012, 2013
 *
 * \file
 * \brief Parses TRANSFAC flat files using boost::spirit qi.
 */


#ifndef BIO_JR_30AUG2011_TRANSFAC_QI_H_
#define BIO_JR_30AUG2011_TRANSFAC_QI_H_



#include "bio/defs.h"

//
// Turn off uninitialised warnings in boost code.
//
#ifdef __GNUC__
# pragma GCC diagnostic push
# pragma GCC diagnostic ignored "-Wuninitialized"
#endif //__GNUC__

#define BOOST_SPIRIT_USE_PHOENIX_V3
#include <boost/spirit/include/phoenix.hpp>
#include <boost/spirit/include/qi.hpp>

#include "bio/spirit/transfac_fusion.h"

#include <boost/spirit/repository/include/qi_kwd.hpp>
#include <boost/spirit/repository/include/qi_keywords.hpp>
#include <boost/fusion/algorithm/transformation/push_back.hpp>
#include <boost/fusion/include/push_back.hpp>
#include <boost/range/algorithm/for_each.hpp>
#include <boost/array.hpp>
#include <boost/range/iterator.hpp>

//
// Turn warnings back on
//
#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif //__GNUC__

#include <iostream>





BIO_NS_START
namespace spirit {

/// Specialised to select correct parser for Entry type.
template< typename Entry, typename Iterator >
struct get_entry_parser;


namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;
namespace phoenix = boost::phoenix;
namespace repository = boost::spirit::repository;

/** Parse a TransData. */
template < typename Iterator >
struct transdata_parser : qi::grammar< Iterator, TransData() >
{

    transdata_parser() : transdata_parser::base_type( start )
    {
        using namespace qi;

        start =
            (lit("CH") | lit("Pathway"))[ _val = PATHWAY_DATA ]
            | lit("Cell")[ _val = CELL_DATA ]
            | (lit("Compel") | lit("C"))[ _val = COMPEL_DATA ]
            | (lit("ev") | lit("Evidence"))[ _val = EVIDENCE_DATA ]
            | (lit("FR") | lit("Fragment"))[ _val = FRAGMENT_DATA ]
            | (lit("Gene") | lit("G"))[ _val = GENE_DATA ]
            | (lit("MO") | lit("Molecule"))[ _val = MOLECULE_DATA ]
            | (lit("Matrix") | lit("M"))[ _val = MATRIX_DATA ]
            | (lit("RE") | lit("Reference"))[ _val = REFERENCE_DATA ]
            | (lit("XN") | lit("Reaction"))[ _val = REACTION_DATA ]
            | (lit("R") | lit("Site"))[ _val = SITE_DATA ]
            | (lit("s") | lit("S"))[ _val = S_DATA ]
            | (lit("T") | lit("Factor"))[ _val = FACTOR_DATA ]
            | eps[ _val = CELL_DATA ]
            ;
    }

    qi::rule< Iterator, TransData() > start;
};



/** Parse a TableLink. */
template < typename Iterator >
struct tablelink_parser : qi::grammar< Iterator, TableLink() >
{
    tablelink_parser() : tablelink_parser::base_type( start, "tablelink" )
    {
        using qi::int_;

        start = transdata >> int_;
    }

    qi::rule< Iterator, TableLink() > start;
    transdata_parser< Iterator > transdata;
};


/** Parse an Identifier. */
template < typename Iterator >
struct identifier_parser : qi::grammar< Iterator, Identifier() >
{
    identifier_parser() : identifier_parser::base_type( start, "identifier" )
    {
        using namespace qi;

        value = *(char_ - '$' - '_' - eol);
        value2 = *(char_ - '$' - eol);
        start =
            hold[(value >> "$" >> -lit("_") >> value >> -("_" >> value2))]
            |
            (attr("") >> value2 >> attr(""))
            ;
    }

    qi::rule< Iterator, Identifier() > start;
    qi::rule< Iterator, std::string() > value;
    qi::rule< Iterator, std::string() > value2;
};


/** Parse a sequence. */
template < typename Iterator >
struct sequence_parser : qi::grammar< Iterator, std::string() >
{
    sequence_parser() : sequence_parser::base_type( start )
    {
        using namespace qi;

        subseq = *(char_ - '.' - eol);
        start = subseq >> *(char_('.') >> subseq) >> omit[-char_('.')];
    }

    qi::rule< Iterator, std::string() > start;
    qi::rule< Iterator, std::string() > subseq;
};


inline
void
print_factorlink( FactorLink & attr ) {
    using namespace boost::phoenix::arg_names;
    std::cout << attr.link << "; " << attr.quality << "; " << attr.name << "; " << attr.sites_included;
    boost::for_each( attr.species, std::cout << "; " << arg1 );
    std::cout << "\n";
}

inline
void
add_species( FactorLink & link, std::string const & value ) {
    //std::cout << "; " << value;
    link.species.push_back( value );
}

inline
void
set_cellular_source( FactorLink & link, std::string const & value ) {
    //std::cout << "; " << value;
    link.cellular_source = value;
}

inline
void
set_quality( FactorLink & link, int value ) {
    //std::cout << "; " << value;
    link.quality = value;
}

inline
void
set_sites_included( FactorLink & link, std::string const & value ) {
    //std::cout << "; " << value;
    if( value == "no." ) {
        link.sites_included = false;
    } else if ( value == "yes." ) {
        link.sites_included = true;
    } else if ( value == "unknown." ) {
        link.sites_included = false;
    } else {
        throw std::logic_error( BIO_MAKE_STRING( "Unknown sites included value: " << value ) );
    }
}


/** Parse a BF entry. */
template < typename Iterator >
struct factorlink_parser : qi::grammar< Iterator, ascii::blank_type, FactorLink() >
{
    factorlink_parser() : factorlink_parser::base_type( start, "factorlink" )
    {
        using namespace qi;
        using namespace keys;
        using qi::_1;
        using boost::phoenix::bind;

        value = lexeme[+(char_ - ";" - eol)]; // a semi-colon or end-of-line terminated value
        start =
            tablelink // the table link
            >> +lit(";")
            >> lexeme[ ( +( char_ - ";" ) ) ] // the name
            >> omit[
                +(
                    lit(";")
                    > ( lit("Cellular source:") >> value[boost::phoenix::bind(set_cellular_source, _val, _1)]
                    |   lit("Species:") >> value[boost::phoenix::bind(add_species, _val, _1)]
                    |   lit("Quality:") >> int_[boost::phoenix::bind(set_quality, _val, _1)]
                    |   lit("site(s) included:") >> value[boost::phoenix::bind(set_sites_included, _val, _1)]
                    )
                )
            ]
            ;
    }

    qi::rule< Iterator, ascii::blank_type, FactorLink() > start;
    qi::rule< Iterator, std::string() > value;
    tablelink_parser< Iterator > tablelink;
};


/**
 * Parse a BS entry.
 *
 * BS  GCTGCAGGTGTTCT; R05108; 3; 15; 17; p.
 */
template < typename Iterator >
struct bindingsite_parser : qi::grammar< Iterator, ascii::blank_type, AlignDesc() >
{
    qi::rule< Iterator, ascii::blank_type, AlignDesc() > start;
    qi::rule< Iterator, ascii::blank_type, char() > base;
    tablelink_parser< Iterator > tablelink;

    bindingsite_parser() : bindingsite_parser::base_type( start, "bindingsite" )
    {
        using namespace qi;

        base =
            (   char_('N')
            |   char_('A')
            |   char_('C')
            |   char_('G')
            |   char_('T')
            |   char_('n')
            |   char_('a')
            |   char_('c')
            |   char_('g')
            |   char_('t')
            |   char_('N')
            |   char_('n')
            |   (omit[alpha] >> attr('N')) // replace anything else with Ns
            )
            ;
        start =
            +base
            > lit(";")
            > tablelink
            > -(lit(",") > tablelink)
            > lit(";")
            > int_ // start
            > lit(";")
            > int_ // length
            > lit(";")
            > -(int_ > *(lit(",") > int_)) // gaps
            > lit(";")
            > ( // orientation
                ('p' >> attr(true))
                | ('n' >> attr(false))
            )
            > '.'
            ;
    }
};


template< typename Iterator >
struct db_ref_string_int_parser : qi::grammar< Iterator, db_ref() >
{
    qi::rule< Iterator, db_ref() > start;
    db_ref_string_int_parser(Database db) : db_ref_string_int_parser::base_type( start )
    {
        using namespace qi;
        start =
            attr(db)
            >> +alpha
            >> int_
            ;
    }
};


template < typename Iterator >
struct affyprobe_parser : qi::grammar< Iterator, db_ref() >
{
    qi::rule< Iterator, db_ref() > start;
    qi::rule< Iterator, std::string() > table;
    qi::rule< Iterator, std::string() > table_end;
    affyprobe_parser() : affyprobe_parser::base_type( start )
    {
        using namespace qi;
        table_end = string("_at") | table;
        table = +(char_ - '_') >> table_end;
        start =
            attr(AFFY_PROBE_DB)
            >> table
            >> attr(0)
            ;
    }
};


template < typename Iterator >
struct dip_parser : qi::grammar< Iterator, db_ref() >
{
    qi::rule< Iterator, db_ref() > start;
    dip_parser() : dip_parser::base_type( start )
    {
        using namespace qi;
        start =
            attr(DIP_DB)
            >> lit("DIP:(") >> +alnum >> lit(")")
            >> attr(0)
            ;
    }
};


template < typename Iterator >
struct embl_parser : qi::grammar< Iterator, db_ref() >
{
    qi::rule< Iterator, std::string() > table;
    qi::rule< Iterator, db_ref() > start;
    embl_parser() : embl_parser::base_type( start )
    {
        using namespace qi;
        table = (
                hold[repeat(4)[upper] >> (repeat(4)[digit]|repeat(8,10)[digit])]
                |
                hold[repeat(3)[upper] >> repeat(5)[digit]]
                |
                hold[repeat(2)[alpha] >> repeat(6)[digit]]
                |
                hold[repeat(1)[upper] >> repeat(5)[digit]]
            )
            >> -(char_('.') >> +digit)
            ;
        start = attr(EMBL_DB) >> table >> attr(int(0));
    }
};


template < typename Iterator >
struct ensembl_parser : db_ref_string_int_parser< Iterator >
{
    ensembl_parser() : db_ref_string_int_parser< Iterator >( ENSEMBL_DB ) { }
};


template < typename Iterator >
struct entrezgene_parser : qi::grammar< Iterator, db_ref() >
{
    qi::rule< Iterator, db_ref() > start;
    entrezgene_parser() : entrezgene_parser::base_type( start )
    {
        using namespace qi;
        start =
            attr(ENTREZ_GENE_DB)
            >> attr("")
            >> int_
            ;
    }
};


template < typename Iterator >
struct entrezprotein_parser : qi::grammar< Iterator, db_ref() >
{
    qi::rule< Iterator, db_ref() > start;
    entrezprotein_parser() : entrezprotein_parser::base_type( start )
    {
        using namespace qi;
        start =
            attr(ENTREZ_PROTEIN_DB)
            >> attr("")
            >> int_
            ;
    }
};


template < typename Iterator >
struct jaspar_parser : db_ref_string_int_parser< Iterator >
{
    jaspar_parser() : db_ref_string_int_parser< Iterator >( JASPAR_DB ) { }
};


template < typename Iterator >
struct flybase_parser : qi::grammar< Iterator, db_ref() >
{
    qi::rule< Iterator, db_ref() > start;
    flybase_parser() : flybase_parser::base_type( start )
    {
        using namespace qi;
        start =
            attr(FLYBASE_DB)
            >> lit("FBgn")
            >> attr("")
            >> int_
            ;
    }
};


template < typename Iterator >
struct prosite_parser : qi::grammar< Iterator, db_ref() >
{
    qi::rule< Iterator, db_ref() > start;
    prosite_parser() : prosite_parser::base_type( start )
    {
        using namespace qi;
        start =
            attr(PROSITE_DB)
            >> lit("PS")
            >> attr("")
            >> int_
            ;
    }
};


template < typename Iterator >
struct refseq_parser : qi::grammar< Iterator, db_ref() >
{
    qi::rule< Iterator, db_ref() > start;
    refseq_parser() : refseq_parser::base_type( start )
    {
        using namespace qi;
        start =
            attr(REFSEQ_DB)
            >> +upper
            >> lit("_")
            >> int_
            ;
    }
};


template < typename Iterator >
struct swissprot_parser : qi::grammar< Iterator, db_ref() >
{
    qi::rule< Iterator, db_ref() > start;
    swissprot_parser() : swissprot_parser::base_type( start )
    {
        using namespace qi;
        start =
            attr(SWISSPROT_DB)
            >> +alnum
            >> attr(0)
            ;
    }
};


template < typename Iterator >
struct transcompel_parser : db_ref_string_int_parser< Iterator >
{
    transcompel_parser() : db_ref_string_int_parser< Iterator >( TRANSCOMPEL_DB ) { }
};


template < typename Iterator >
struct transfac_parser : db_ref_string_int_parser< Iterator >
{
    transfac_parser() : db_ref_string_int_parser< Iterator >( TRANSFAC_DB ) { }
};


template < typename Iterator >
struct transpath_parser : db_ref_string_int_parser< Iterator >
{
    transpath_parser() : db_ref_string_int_parser< Iterator >( TRANSPATH_DB ) { }
};


template < typename Iterator >
struct transpro_parser : qi::grammar< Iterator, db_ref() >
{
    qi::rule< Iterator, db_ref() > start;
    transpro_parser() : transpro_parser::base_type( start )
    {
        using namespace qi;
        start =
            attr(TRANSPRO_DB)
            >> +(char_ - '.')
            >> attr(0)
            ;
    }
};


/** Parse a database reference. */
template < typename Iterator >
struct dbref_parser : qi::grammar< Iterator, db_ref(), ascii::blank_type >
{
    qi::rule< Iterator, db_ref(), ascii::blank_type > start;
    qi::rule< Iterator, db_ref() > unknown; // an unknown dbref
    qi::rule< Iterator > ignore; // ignore until eol
    qi::rule< Iterator > notparseddb; // databases that won't be parsed
    qi::rule< Iterator, db_ref(), ascii::blank_type >
        affyline, dipline, emblline, ensemblline, entrezgeneline,
        entrezproteinline, epdline, flybaseline, jasparline, prositeline, refseqline,
        swissprotline, transcompelline, transfacline, transpathline, transproline,
        uniprotline, notparsedline
        ;
    affyprobe_parser< Iterator > affyprobe;
    dip_parser< Iterator > dip;
    ensembl_parser< Iterator > ensembl;
    embl_parser< Iterator > embl;
    entrezgene_parser< Iterator > entrezgene;
    entrezprotein_parser< Iterator > entrezprotein;
    flybase_parser< Iterator > flybase;
    jaspar_parser< Iterator > jaspar;
    prosite_parser< Iterator > prosite;
    refseq_parser< Iterator > refseq;
    swissprot_parser< Iterator > swissprot;
    transcompel_parser< Iterator > transcompel;
    transfac_parser< Iterator > transfac;
    transpath_parser< Iterator > transpath;
    transpro_parser< Iterator > transpro;

    dbref_parser() : dbref_parser::base_type( start, "dbref" )
    {
        using namespace qi;

        ignore = omit[*(char_ - eol)];
        unknown = ignore >> attr(db_ref());
        notparseddb =
                lit("BKL:")
            |   lit("BRENDA:")
            |   lit("CLDB:")
            |   lit("DATF:")
            |   lit("EPD:")
            |    lit("HGNC:")
            |   lit("MGI:")
            |   lit("MIRBASE:")
            |   lit("OMIM:")
            |   lit("PATHODB:")
            |   lit("PDB:")
            |   lit("PIR:")
            |   lit("RGD:")
            |   lit("RSNP:")
            |   lit("SMARTDB:")
            |   lit("TRRD:")
            |   lit("UNIGENE:")
            |   lit("WB:")
            |   lit("ZFIN:")
            ;
        affyline = lit("AFFYMETRIX:") >> (-lit("Affymetrix:")) > unknown;
        //            (lit("AFFYMETRIX:") >> (-lit("Affymetrix:")) > affyprobe > ignore) // affy parser doesn't work yet
        dipline = lit("DIP:") > dip > ignore;
        emblline = lit("EMBL:") > -lit( ';' ) > (embl|refseq|swissprot) > ignore; // some other entries are given as EMBL
        ensemblline = lit("ENSEMBL:") > ensembl > ignore;
        entrezgeneline = lit("ENTREZGENE:") > entrezgene > ignore;
        entrezproteinline = lit("ENTREZPROTEIN:") > entrezprotein > ignore;
        flybaseline = (lit("FLYBASE:")|lit("FB:")|lit("Flybase:")) > flybase > ignore;
        jasparline = lit("JASPAR:") > jaspar > ignore;
        prositeline = lit("PROSITE:") > prosite > ignore;
        refseqline = lit("REFSEQ:") > refseq > ignore;
        swissprotline = lit("SWISSPROT:") > swissprot > ignore;
        transcompelline = lit("TRANSCOMPEL:") > transcompel > ignore;
        transfacline = lit("TRANSFAC:") > transfac > ignore;
        transpathline = lit("TRANSPATH:") > transpath > ignore;
        transproline = lit("TRANSPRO:") > transpro > ignore;
        uniprotline = lit("UNIPROT:") > swissprot > ignore;
        notparsedline = notparseddb > unknown;
        start =
                affyline
            |   dipline
            |   emblline
            |   ensemblline
            |   entrezgeneline
            |   entrezproteinline
            |   flybaseline
            |   jasparline
            |   prositeline
            |   refseqline
            |   swissprotline
            |   transcompelline
            |   transfacline
            |   transpathline
            |   transproline
            |   uniprotline
            |   notparsedline
            |   lit("DOES NOT MATCH") // comment this line out and uncomment the following line
            //unknown
            ;
    }
};


/** Parse a database. */
template < typename Iterator >
struct db_parser : qi::grammar< Iterator, Database() >
{
    db_parser() : db_parser::base_type( start, "db" )
    {
        using namespace qi;

        start =
            lit("AFFYMETRIX")[ _val = AFFY_PROBE_DB ]
            | lit("BKL")[ _val = BKL_DB ]
            | lit("DIP")[ _val = DIP_DB ]
            | lit("EMBL")[ _val = EMBL_DB ]
            | lit("ENSEMBL")[ _val = ENSEMBL_DB ]
            | lit("ENTREZGENE")[ _val = ENTREZ_GENE_DB ]
            | lit("ENTREZPROTEIN")[ _val = ENTREZ_PROTEIN_DB ]
            | lit("EPD")[ _val = EPD_DB ]
            | (lit("FLYBASE")|lit("FB"))[ _val = FLYBASE_DB ]
            | lit("JASPAR")[ _val = JASPAR_DB ]
            | lit("MGI")[ _val = MGI_DB ]
            | lit("PATHODB")[ _val = PATHO_DB ]
            | lit("PDB")[ _val = PDB_DB ]
            | lit("PIR")[ _val = PIR_DB ]
            | lit("REFSEQ")[ _val = REFSEQ_DB ]
            | lit("RGD")[ _val = RGD_DB ]
            | lit("RSNP")[ _val = RSNP_DB ]
            | lit("SMARTDB")[ _val = SMART_DB ]
            | (lit("SWISSPROT")|lit("UNIPROT"))[ _val = SWISSPROT_DB ]
            | lit("TRANSCOMPEL")[ _val = TRANSCOMPEL_DB ]
            | lit("TRANSFAC")[ _val = TRANSFAC_DB ]
            | lit("TRANSPATH")[ _val = TRANSPATH_DB ]
            | lit("TRANSPRO")[ _val = TRANSPRO_DB ]
            | lit("UNIGENE")[ _val = UNIGENE_DB ]
            | lit("WB")[ _val = WORMBASE_DB ]
            | lit("ZFIN")[ _val = ZFIN_DB ]
            ;
    }

    qi::rule< Iterator, Database() > start;
};


inline
void
parse_dbref( db_ref & value, Database db, std::string const & acc ) {
    value = parse_db_ref_as( acc, db );
}

typedef boost::array< BIO_NS::float_t, 4 > pssm_counts_type;

/** Parse a PssmEntry. */
template < typename Iterator >
struct pssmentry_parser : qi::grammar< Iterator, ascii::blank_type, PssmEntry(), qi::locals< pssm_counts_type > >
{
    qi::rule< Iterator, ascii::blank_type, PssmEntry(), qi::locals< pssm_counts_type > > start;

    pssmentry_parser() : pssmentry_parser::base_type( start, "pssmentry" )
    {
        using namespace qi;
        using namespace phoenix;
        using qi::_1;

        start =
            double_[ _a[ 0 ] = _1 ]
            > double_[ _a[ 1 ] = _1 ]
            > double_[ _a[ 2 ] = _1 ]
            > double_[ _a[ 3 ] = _1 ]
            > eps[ _val = construct< PssmEntry >( phoenix::bind( &pssm_counts_type::elems, _a ) ) ]
            ;
    }
};


/** Parse a Pssm section. */
template < typename Iterator >
struct pssm_section_parser : qi::grammar< Iterator, ascii::blank_type, PssmAndConsensus() >
{
    pssmentry_parser< Iterator > pssmentry;
    qi::rule< Iterator, ascii::blank_type, PssmAndConsensus() > start;

    pssm_section_parser() : pssm_section_parser::base_type( start, "pssmsection" )
    {
        using namespace qi;
        using namespace phoenix;
        using boost::phoenix::push_back;
        using boost::phoenix::bind;
        using qi::_1;

        start =
            lit("P0      A      C      G      T") > eol
            > *(  omit[ repeat(2)[digit] ]
                > pssmentry[ push_back( bind( &PssmAndConsensus::pssm, _val ), _1 ) ]
                > upper[ push_back( bind( &PssmAndConsensus::consensus, _val ), _1 ) ]
                > eol
            )
            ;
    }
};


/**
 * Parses a line. Has same attribute type as parser, P
 */
template< typename Iterator, typename P, bool ReturnTag=false >
struct biobaseline_parser
: qi::grammar<
      Iterator,
      typename boost::mpl::deref< typename boost::mpl::begin< typename P::sig_type >::type >::type,
      ascii::blank_type
>
{
    qi::rule<
        Iterator,
        typename boost::mpl::deref< typename boost::mpl::begin< typename P::sig_type >::type >::type,
        ascii::blank_type
    > start;
    P p;
    biobaseline_parser( char const * tag ) : biobaseline_parser::base_type( start ) {
        using namespace qi;
        if( ReturnTag ) {
            start = string(tag) >> p;
        } else {
            start = lit(tag) >> p;
        }
    }
};


template< typename Range >
boost::optional< int >
parse_matrix_number_of_sites( Range const & basis ) {
    boost::optional< int > result;
    int i;
    if( qi::parse(
            boost::begin( basis ),
            boost::end( basis ),
            qi::int_,
            i
        )
    ) {
        result = i;
    }
    return result;
}


/** Tidy up any loose ends. */
inline
void
finalise_matrix( Matrix & matrix ) {

    //
    // set up the number of sites
    //
    // try to parse number of sites from basis string
    const boost::optional< int > num_sites = parse_matrix_number_of_sites( matrix.matrix_basis );

    // if we got a number of sites then use it
    if( num_sites ) {

        matrix.number_of_sites = *num_sites;

    } else {

        // otherwise see how many observations are in PSSM
        const float_t num_obs = matrix.pssm.empty() ? 0. : matrix.pssm[0].get_num_observations();

        // if we don't have many observations in our PSSM counts perhaps they are frequencies
        if( num_obs < 2. ) {
            // so we use a default number of sites
            matrix.number_of_sites = ANALOGUE_SITE_EQUIVALENT;
        } else {
            // otherwise we take the number of observations as the number of sites
            matrix.number_of_sites = num_obs;
        }
    }
}



/**
 * Base class for parsers of biobase table entries. Uses CRTP idiom.
 */
template< typename Iterator >
struct entry_parser
{
    // grammars
    tablelink_parser< Iterator > tablelink;
    identifier_parser< Iterator > identifier;
    factorlink_parser< Iterator > factorlink;
    dbref_parser< Iterator > dbref;
    bindingsite_parser< Iterator > bindingsite;
    sequence_parser< Iterator > sequence;

    // rules
    qi::rule< Iterator, std::string() > rest_of_line, name_section, description_section,
        concatenated_sequence_section, ///< Concatenates all sequences together
        first_sequence_section, ///< Just returns first sequence
        multiline_sequence,
        sequence_line;
    qi::rule< Iterator > section_separator, ignored_rest_of_line, ignored_section;
    qi::rule< Iterator, std::set< std::string >() > comma_separated_string_set;
    qi::rule< Iterator, std::vector< std::string >() > comma_separated_string_vec,
        individual_sequences_section; ///< Returns vector of sequences
    qi::rule< Iterator, TableLink( char const * ) > tablelink_line;
    qi::rule< Iterator, std::string( char const * ) > string_line;
    qi::rule< Iterator, TableLink() > accession_section;
    qi::rule< Iterator, DatabaseRefVec() > dbref_section;
    qi::rule< Iterator, AlignDescPtr() > binding_site_ptr;
    qi::rule< Iterator, AlignDescList() > binding_site_section;
    qi::rule< Iterator, FactorLinkPtr() > factor_link_ptr;
    qi::rule< Iterator, FactorLinkList() > factor_link_section;
    qi::rule< Iterator, Identifier() > identifier_section;

    // symbols
    qi::symbols< char > ignored_key; ///< Ignored keys

    entry_parser()
    {
        using namespace qi;
        using namespace phoenix;
        using qi::_1;

        //
        // Helper rules
        //
        rest_of_line                   = *( char_ - eol );
        ignored_rest_of_line           = omit[ rest_of_line > eol ];
        section_separator              = lit("XX") > eol;
        binding_site_ptr               = skip( ascii::blank )[
                                           bindingsite[ _val = construct< AlignDescPtr >( new_< AlignDesc >( _1 ) ) ]
                                       ];
        factor_link_ptr                = skip( ascii::blank )[
                                           factorlink[ _val = construct< FactorLinkPtr >( new_< FactorLink >( _1 ) ) ]
                                       ];
        comma_separated_string_set     = *( char_ - ';' - eol ) % ( lit( ';' ) > -lit( ' ' ) );
        comma_separated_string_vec     = *( char_ - ';' - eol ) % ( lit( ';' ) > -lit( ' ' ) );
        tablelink_line                 = lit( _r1 ) > tablelink > ignored_rest_of_line;
        string_line                    = lit( _r1 ) > rest_of_line > eol;
        sequence_line                  =
                lit( "SQ  " )
            >   *( char_ - ( char_( '.' ) >> eol ) - eol )
            >   -omit[ char_( '.' ) ]
            >   eol;
        multiline_sequence             = ( lit( "SQ  " ) > *( char_ - eol - ('.' >> eol) ) ) % eol > '.' > eol;


        //
        // Sections
        //
        ignored_section                = +( ignored_key > ignored_rest_of_line );
        accession_section              = tablelink_line( val( "AC  " ) );
        identifier_section             = lit( "ID  " ) > identifier > eol;
        name_section                   = string_line( val( "NA  " ) );
        description_section            = string_line( val( "DE  " ) );
        concatenated_sequence_section  = +sequence_line;
        first_sequence_section         = sequence_line > omit[ *sequence_line ];
        individual_sequences_section   = +sequence_line;
        binding_site_section           = +( lit( "BS  " ) > binding_site_ptr > eol );
        dbref_section                  = +( lit( "DR  " ) > skip( ascii::blank )[ dbref ] > eol );
        factor_link_section            = +( lit( "BF  " ) > factor_link_ptr > eol );

        //
        // Name the rules for debugging
        //
        rest_of_line.name( "rest of line" );
        section_separator.name( "section separator" );
        ignored_rest_of_line.name( "ignored rest of line" );
        accession_section.name( "accession_section" );
        dbref_section.name( "dbref_section" );
        binding_site_ptr.name( "binding_site_ptr" );
        binding_site_section.name( "binding_site_section" );
        factor_link_ptr.name( "factor_link_ptr" );
        factor_link_section.name( "factor_link_section" );
        name_section.name( "name_section" );
        description_section.name( "description_section" );
        concatenated_sequence_section.name( "concatenated_sequence_section" );
        individual_sequences_section.name( "individual_sequences_section" );

//        debug( binding_site_ptr );
//        debug( binding_site_section );
    }
};




/**
 * Parses an entry in the pathway table.
 */
template< typename Iterator >
struct pathway_parser
: qi::grammar< Iterator, Pathway::ptr_t() >
, entry_parser< Iterator >
{
    typedef Pathway parsed_type;
    typedef typename parsed_type::ptr_t parsed_ptr;

    // rules
    qi::rule< Iterator, parsed_ptr() > start;

    pathway_parser() : pathway_parser::base_type( start ) {
        using namespace phoenix;
        using namespace qi;
        using namespace keys;
        using qi::_1;

        // The keys that are used are commented out.
        this->ignored_key.add
//            ( "AC  " ) // Accession number
            ( "DT  " ) // Created/Updated
            ( "CO  " ) // Copyright information
//            ( "CE  " ) // Composite element
            ( "EX  " ) // Experiment type
            ( "CN  " ) // Experiment conclusion
            ( "CT  " ) // Cell line
            ( "CC  " ) // Evidence comment
            ( "IN  " ) // Interaction
            ( "FA  " ) // Factor
            ( "FD  " ) // DNA-binding domain
            ( "FF  " ) // Effect on transcription
            ( "FC  " ) // Factor comment
//            ( "DR  " ) // Link to external databases: TRANSFAC
            ( "RN  " ) // Reference number
            ( "RX  " ) // Reference ...
            ( "RA  " ) // Author(s)
            ( "RT  " ) // Title
            ( "RL  " ) // Journal
            ;

        start = // parse each section
                eps[ _val = construct< parsed_ptr >( new_< parsed_type >() ) ] // construct the shared_ptr
            >>  +(  (   this->ignored_section // ignore
                    |   this->accession_section             [ phoenix::bind( &parsed_type::accession_number , *_val ) = _1 ]
                    ) > -this->section_separator // sometimes Biobase forget the separator so make it optional
                )
            ;
    }
};

/// Partial specialisation to select correct parser
template< typename Iterator >
struct get_entry_parser< Pathway, Iterator > {
    typedef pathway_parser< Iterator > type;
};



/**
 * Parses an entry in the moleule table.
 */
template< typename Iterator >
struct molecule_parser
: qi::grammar< Iterator, Molecule::ptr_t() >
, entry_parser< Iterator >
{
    typedef Molecule parsed_type;
    typedef typename parsed_type::ptr_t parsed_ptr;

    // rules
    qi::rule< Iterator, parsed_ptr() > start;

    molecule_parser() : molecule_parser::base_type( start ) {
        using namespace phoenix;
        using namespace qi;
        using namespace keys;
        using qi::_1;

        // The keys that are used are commented out.
        this->ignored_key.add
//            ( "AC  " ) // Accession number
            ( "DT  " ) // Created/Updated
            ( "CO  " ) // Copyright information
//            ( "CE  " ) // Composite element
            ( "EX  " ) // Experiment type
            ( "CN  " ) // Experiment conclusion
            ( "CT  " ) // Cell line
            ( "CC  " ) // Evidence comment
            ( "IN  " ) // Interaction
            ( "FA  " ) // Factor
            ( "FD  " ) // DNA-binding domain
            ( "FF  " ) // Effect on transcription
            ( "FC  " ) // Factor comment
//            ( "DR  " ) // Link to external databases: TRANSFAC
            ( "RN  " ) // Reference number
            ( "RX  " ) // Reference ...
            ( "RA  " ) // Author(s)
            ( "RT  " ) // Title
            ( "RL  " ) // Journal
            ;

        start = // parse each section
                eps[ _val = construct< parsed_ptr >( new_< parsed_type >() ) ] // construct the shared_ptr
            >>  +(  (   this->ignored_section // ignore
                    |   this->accession_section             [ phoenix::bind( &parsed_type::accession_number , *_val ) = _1 ]
                    ) > -this->section_separator // sometimes Biobase forget the separator so make it optional
                )
            ;
    }
};

/// Partial specialisation to select correct parser
template< typename Iterator >
struct get_entry_parser< Molecule, Iterator > {
    typedef molecule_parser< Iterator > type;
};



/**
 * Parses an entry in the compel table.
 */
template< typename Iterator >
struct compel_parser
: qi::grammar< Iterator, Compel::ptr_t() >
, entry_parser< Iterator >
{
    typedef Compel parsed_type;
    typedef typename parsed_type::ptr_t parsed_ptr;
    typedef std::pair< int, int > from_to_range;

    // rules
    qi::rule< Iterator, parsed_ptr() > start;
    qi::rule< Iterator, ascii::blank_type, from_to_range() > from_to;
    qi::rule< Iterator, ascii::blank_type, CompelBindingSite() > bindingsite;

    qi::symbols< char, CompelType > type;

    compel_parser() : compel_parser::base_type( start ) {
        using namespace phoenix;
        using namespace qi;
        using namespace keys;
        using qi::_1;

        type.add
            ( "synergism",  COMPEL_SYNERGISM    )   //synergism
            ( "antagonism", COMPEL_ANTAGONISM   )   //antagonism
            ;

        // The keys that are used are commented out.
        this->ignored_key.add
//            ( "AC  " ) // Accession number
//            ( "ID  " ) // Identifier
            ( "DT  " ) // Created/Updated
            ( "CO  " ) // Copyright information
//            ( "GE  " ) // TRANSFAC gene
//            ( "SQ  " ) // Sequence
            ( "ST  " ) // Reference point for sequence start
//            ( "PS  " ) // Position
//            ( "DR  " ) // Link to external databases: EMBL
//            ( "NO  " ) // Number of sites
//            ( "BS  " ) // Binding site
//            ( "TY  " ) // Type (synergism/antagonism)
            ( "CL  " ) // Functional classification
            ( "CG  " ) // Compel group
            ( "CM  " ) // Molecular mechanism
//            ( "CC  " ) // Comment
//            ( "EV  " ) // Evidence
            ( "EX  " ) // Experiment type
            ( "CN  " ) // Experiment conclusion
            ( "CT  " ) // Cell line
            ( "RN  " ) // Reference number
            ( "RX  " ) // Reference ...
            ( "RA  " ) // Author(s)
            ( "RT  " ) // Title
            ( "RL  " ) // Journal
            ( "C4  " ) //    C4 zinc finger type of nuclear receptors: GR, PR, ER, RAR, T3R,
            ( "WH  " ) //    HTH, winged helix/fork head (HNF-3, E2F, FAST, RFX families)
            ;

        from_to =
                int_[ phoenix::bind( &from_to_range::first , _val ) = _1 ]
            >   "to"
            >   int_[ phoenix::bind( &from_to_range::second, _val ) = _1 ]
            ;

        bindingsite =
                from_to[
                    phoenix::bind( &CompelBindingSite::start, _val ) = phoenix::bind( &from_to_range::first , _1 ),
                    phoenix::bind( &CompelBindingSite::end  , _val ) = phoenix::bind( &from_to_range::second, _1 )
                ]
            >   ';'
            >   as_string[ +( char_ - ';' ) ][ phoenix::bind( &CompelBindingSite::factor , _val ) = _1 ]
            >   ';'
            >   this->tablelink[ phoenix::bind( &CompelBindingSite::site_link , _val ) = _1 ]
            >   '.'
            ;

        start = // parse each section
                eps[ _val = construct< parsed_ptr >( new_< parsed_type >() ) ] // construct the shared_ptr
            >>  +(  (   this->ignored_section // ignore
                    |   this->accession_section              [ phoenix::bind( &parsed_type::accession_number , *_val ) = _1 ]
                    |   this->string_line( val( "ID  " ) )   //ignore the identifier
                    |   this->dbref_section                  [ phoenix::bind( &parsed_type::database_refs    , *_val ) = _1 ]
                    |   this->concatenated_sequence_section  [ phoenix::bind( &parsed_type::sequence         , *_val ) = _1 ]
                    |   this->dbref_section                  [ phoenix::bind( &parsed_type::database_refs    , *_val ) = _1 ]
                    |   +this->string_line( val( "CC  " ) )  [ phoenix::bind( &parsed_type::comment, *_val ) +=_1 ]
                    |   ( lit( "TY  " ) > type                [ phoenix::bind( &parsed_type::type             , *_val ) = _1 ] > eol )
                    |   (   lit( "NO  " ) > int_ > eol
                        >   +(  lit( "BS  " )
                            >   skip( ascii::blank )[
                                    bindingsite[
                                        phoenix::push_back( phoenix::bind( &parsed_type::binding_sites, *_val ), _1 )
                                    ]
                                ]
                            >   eol
                            )
                        )
                    |   (   this->tablelink_line( val( "GE  " ) )[ phoenix::bind( &parsed_type::gene, *_val ) = _1 ]
                        >   omit[ *this->string_line( val( "GE  " ) ) ]
                        )
                    |   // position section
                        (   lit( "PS  " )
                        >   skip( ascii::blank )[
                                from_to[
                                    phoenix::bind( &parsed_type::begin, *_val ) = phoenix::bind( &from_to_range::first , _1 ),
                                    phoenix::bind( &parsed_type::end  , *_val ) = phoenix::bind( &from_to_range::second, _1 )
                                ]
                            ]
                        >   eol
                        )
                    |   (   lit( "EV  " ) > this->tablelink[ phoenix::push_back( phoenix::bind( &parsed_type::evidences, *_val ), _1 ) ] > eol
                        >   -( lit( "EX  " ) > omit[ +( char_ - eol ) ] > eol )
                        >   -( lit( "CN  " ) > omit[ +( char_ - eol ) ] > eol )
                        )
                    ) > -this->section_separator // sometimes Biobase forget the separator so make it optional
                )
            ;
    }
};

/// Partial specialisation to select correct parser
template< typename Iterator >
struct get_entry_parser< Compel, Iterator > {
    typedef compel_parser< Iterator > type;
};



/**
 * Parses an entry in the evidence table.
 */
template< typename Iterator >
struct evidence_parser
: qi::grammar< Iterator, Evidence::ptr_t() >
, entry_parser< Iterator >
{
    typedef Evidence parsed_type;
    typedef typename parsed_type::ptr_t parsed_ptr;

    // rules
    qi::rule< Iterator, parsed_ptr() > start;

    evidence_parser() : evidence_parser::base_type( start ) {
        using namespace phoenix;
        using namespace qi;
        using namespace keys;
        using qi::_1;

        // The keys that are used are commented out.
        this->ignored_key.add
//            ( "AC  " ) // Accession number
            ( "DT  " ) // Created/Updated
            ( "CO  " ) // Copyright information
//            ( "CE  " ) // Composite element
            ( "EX  " ) // Experiment type
            ( "CN  " ) // Experiment conclusion
            ( "CT  " ) // Cell line
            ( "CC  " ) // Evidence comment
            ( "IN  " ) // Interaction
            ( "FA  " ) // Factor
            ( "FD  " ) // DNA-binding domain
            ( "FF  " ) // Effect on transcription
            ( "FC  " ) // Factor comment
//            ( "DR  " ) // Link to external databases: TRANSFAC
            ( "RN  " ) // Reference number
            ( "RX  " ) // Reference ...
            ( "RA  " ) // Author(s)
            ( "RT  " ) // Title
            ( "RL  " ) // Journal
            ;

        start = // parse each section
                eps[ _val = construct< parsed_ptr >( new_< parsed_type >() ) ] // construct the shared_ptr
            >>  +(  (   this->ignored_section // ignore
                    |   this->accession_section             [ phoenix::bind( &parsed_type::accession_number , *_val ) = _1 ]
                    |   this->tablelink_line( val( "CE  ") )[ phoenix::bind( &parsed_type::composite_element, *_val ) = _1 ]
                    |   this->dbref_section                 [ phoenix::bind( &parsed_type::database_refs    , *_val ) = _1 ]
                    ) > -this->section_separator // sometimes Biobase forget the separator so make it optional
                )
            ;
    }
};

/// Partial specialisation to select correct parser
template< typename Iterator >
struct get_entry_parser< Evidence, Iterator > {
    typedef evidence_parser< Iterator > type;
};



/**
 * Parses an entry in the gene table.
 */
template< typename Iterator >
struct gene_parser
: qi::grammar< Iterator, Gene::ptr_t() >
, entry_parser< Iterator >
{
    typedef Gene parsed_type;
    typedef typename parsed_type::ptr_t parsed_ptr;

    // rules
    qi::rule< Iterator, parsed_ptr() > start;

    gene_parser() : gene_parser::base_type( start ) {
        using namespace phoenix;
        using namespace qi;
        using namespace keys;
        using qi::_1;

        // The keys that are used are commented out.
        this->ignored_key.add
//            ( "AC  " ) // Accession number
            ( "AS  " ) // Accession numbers, secondary
//            ( "ID  " ) // Identifier
            ( "DT  " ) // Created/Updated
            ( "CO  " ) // Copyright
            ( "SD  " ) // Short description
            ( "DE  " ) // Description
            ( "SY  " ) // Synonyms
            ( "OS  " ) // Species
            ( "OC  " ) // Taxonomic classification
            ( "CH  " ) // Chromosomal location
            ( "HG  " ) // Host gene
            ( "IG  " ) // Intronic gene
            ( "BC  " ) // EPD Promoter classification
            ( "RG  " ) // Regulation
            ( "BS  " ) // Binding sites / Binding factors
            ( "BR  " ) // Binding region (ChIP-chip/-Seq)
            ( "CE  " ) // Composite element
            ( "FA  " ) // Encoded factor
//            ( "DR  " ) // External database links
            ( "RN  " ) // Reference number
            ( "RX  " ) // PUBMED; link to PubMed entry.
            ( "RA  " ) // Reference authors
            ( "RT  " ) // Reference title
            ( "RL  " ) // Reference source
            ;

        start = // parse each section
                eps[ _val = construct< parsed_ptr >( new_< parsed_type >() ) ] // construct the shared_ptr
            >>  +(  (   this->ignored_section // ignore
                    |   this->accession_section[ phoenix::bind( &parsed_type::accession_number, *_val ) = _1 ]
                    |   this->dbref_section    [ phoenix::bind( &parsed_type::database_refs, *_val ) = _1 ]
                    |   // name & species in identifier section
                        (   lit( "ID  " )
                        >   as_string[ ( *( char_ - '$' ) ) ][ phoenix::bind( &parsed_type::species, *_val ) = _1 ]
                        >   '$'
                        >   as_string[ ( *( char_ - eol ) ) ][ phoenix::bind( &parsed_type::name   , *_val ) = _1 ]
                        >   eol
                        )
                    ) > -this->section_separator // sometimes Biobase forget the separator so make it optional
                )
            ;
    }
};

/// Partial specialisation to select correct parser
template< typename Iterator >
struct get_entry_parser< Gene, Iterator > {
    typedef gene_parser< Iterator > type;
};



/**
 * Parses a factor entry in a table.
 */
template< typename Iterator >
struct factor_parser
: qi::grammar< Iterator, Factor::ptr_t() >
, entry_parser< Iterator >
{
    typedef Factor parsed_type;
    typedef typename parsed_type::ptr_t parsed_ptr;

    qi::rule< Iterator, parsed_ptr() > start;
    qi::symbols< char, Factor::type > type;

    factor_parser() : factor_parser::base_type( start, "factor" ) {
        using namespace phoenix;
        using namespace qi;
        using namespace keys;
        using qi::_1;

        this->ignored_key.add
            ( "AS  " ) // secondary accession
            ( "DT  " ) // date
            ( "CO  " ) // copyright
            ( "RN  " ) // reference
            ( "RX  " ) // ...
            ( "RA  " ) // ...
            ( "RT  " ) // ...
            ( "RL  " ) // ...
            ( "ID  " ) // identifier
            ( "FT  " ) // feature table
            ( "FF  " ) // functional property
            ( "SF  " ) // structural feature
            ( "BS  " ) // site
            ( "HO  " ) // homolog
            ( "CL  " ) // class
            ( "CP  " ) // Cell specificity (positive)
            ( "CN  " ) // Cell specificity (negative)
            ( "EX  " ) // Expression pattern
            ( "IN  " ) // Interacting factors
            ( "BS  " ) // Binding sites / Regulated genes
            ( "BR  " ) // Binding region (ChIP-chip/-Seq)
            ( "TY  " ) // Type
            ( "SZ  " ) // Size
            ( "SQ  " ) // Sequence
            ( "SC  " ) // Sequence source
            ;

        type.add
            ( "family"         , Factor::family_type )
            ( "family_mod"     , Factor::family_mod_type )
            ( "isogroup"       , Factor::isogroup_type )
            ( "isogroup_mod"   , Factor::isogroup_mod_type )
            ( "complex"        , Factor::complex_type )
            ( "complex_mod"    , Factor::complex_mod_type )
            ( "basic"          , Factor::basic_type )
            ( "basic_mod"      , Factor::basic_mod_type )
            ( "miRNA basic"    , Factor::miRNA_type )
            ( "pre-miRNA basic", Factor::pre_miRNA_type )
            ;

        start = // all the sections
                eps[ _val = construct< parsed_ptr >( new_< parsed_type >() ) ]
            >>  +(  (   this->ignored_section // ignore
                    |   this->accession_section[ phoenix::bind( &parsed_type::accession_number, *_val ) = _1 ]
                    |   this->dbref_section[ phoenix::bind( &parsed_type::database_refs, *_val ) = _1 ]
                    |   ( lit( "TY  " ) > type[ phoenix::bind( &parsed_type::_type, *_val ) = _1 ] > '.' > eol )
                    |   (   lit( "SY  " )
                        >   this->comma_separated_string_set[ phoenix::bind( &parsed_type::synonyms, *_val ) = _1 ]
                        >   eol
                        )
                    |   (   lit( "OS  " )
                        >   this->ignored_rest_of_line
                        >   -(  lit( "OC  " )
                            >   this->comma_separated_string_vec[ phoenix::bind( &parsed_type::taxonomies, *_val ) = _1 ]
                            >   eol
                            )
                        )
                    |   this->tablelink_line( val( "GE  " ) )[ phoenix::bind( &parsed_type::gene, *_val ) = _1 ]
                    |   this->string_line( val( "FA  " ) )[ phoenix::bind( &parsed_type::name, *_val ) = _1 ]
                    |   ( +this->tablelink_line( val( "MX  " ) ) )[ phoenix::bind( &parsed_type::matrices, *_val ) = _1 ]
                    |   ( +this->tablelink_line( val( "ST  " ) ) )[ phoenix::bind( &parsed_type::subunits, *_val ) = _1 ]
                    |   ( +this->tablelink_line( val( "CX  " ) ) )[ phoenix::bind( &parsed_type::complexes, *_val ) = _1 ]
                    |   ( +this->tablelink_line( val( "HC  " ) ) )[ phoenix::bind( &parsed_type::sub_families, *_val ) = _1 ]
                    |   ( +this->tablelink_line( val( "HP  " ) ) )[ phoenix::bind( &parsed_type::super_families, *_val ) = _1 ]
                    ) > -this->section_separator // sometimes Biobase forget the separator so make it optional
                )
            ;
//        debug( start );
//        debug( this->accession_section );
//        debug( this->identifier_section );
//        debug( this->binding_site_section );
//        debug( this->dbref_section );
//        debug( this->section_separator );
    }
};

/// Partial specialisation to select correct parser
template< typename Iterator >
struct get_entry_parser< Factor, Iterator > {
    typedef factor_parser< Iterator > type;
};



/**
 * Parses a matrix entry in a table.
 */
template< typename Iterator >
struct matrix_parser
: qi::grammar< Iterator, Matrix::ptr_t() >
, entry_parser< Iterator >
{
    typedef Matrix parsed_type;
    typedef typename parsed_type::ptr_t parsed_ptr;

    // grammars
    pssm_section_parser< Iterator > pssmsection;

    // rules
    qi::rule< Iterator, parsed_ptr() > start;

    matrix_parser() : matrix_parser::base_type( start ) {
        using namespace phoenix;
        using namespace qi;
        using namespace keys;
        using qi::_1;

        // The keys that are used are commented out.
        this->ignored_key.add
//            ( "AC  " ) // Accession number
            ( "AS  " ) // Accession numbers, secondary
//            ( "ID  " ) // Identifier
            ( "DT  " ) // Created/Updated
            ( "CO  " ) // Copyright
//            ( "NA  " ) // Name
//            ( "DE  " ) // Factor description
            ( "TY  " ) // Type
            ( "OS  " ) // Species/Taxon
            ( "OC  " ) // Taxonomic classification
            ( "HP  " ) // Superfamilies
            ( "HC  " ) // Subfamilies
//            ( "BF  " ) // Binding Factors
//            ( "P0  " ) // Binding Matrix
//            ( "BA  " ) // Basis
//            ( "BS  " ) // Binding sites
            ( "CC  " ) // Comments
            ( "DR  " ) // External database links
            ( "OV  " ) // Older version
            ( "PV  " ) // Preferred version
            ( "RN  " ) // Reference number
            ( "RX  " ) // PUBMED; link to PubMed entry.
            ( "RA  " ) // Reference authors
            ( "RT  " ) // Reference title
            ( "RL  " ) // Reference source
            ;

        start = // all the sections
                eps[ _val = construct< parsed_ptr >( new_< parsed_type >() ) ]
            >>  +(  (   this->ignored_section // ignore
                    |   this->accession_section[ phoenix::bind( &parsed_type::accession_number, *_val ) = _1 ]
                    |   this->identifier_section[ phoenix::bind( &parsed_type::id, *_val ) = _1 ]
                    |   this->name_section[ phoenix::bind( &parsed_type::factor_name, *_val ) = _1 ]
                    |   this->description_section[ phoenix::bind( &parsed_type::description, *_val ) = _1 ]
                    |   this->factor_link_section[ phoenix::bind( &parsed_type::factor_links, *_val ) = _1 ]
                    |   this->binding_site_section[ phoenix::bind( &parsed_type::align_descs, *_val ) = _1 ]
                    |   this->string_line( val( "BA  " ) )[ phoenix::bind( &parsed_type::matrix_basis, *_val ) = _1 ]
                    |   skip( ascii::blank )[
                            pssmsection
                                [   phoenix::bind( &parsed_type::pssm, *_val ) = phoenix::bind( &PssmAndConsensus::pssm, _1 ),
                                    phoenix::bind( &parsed_type::consensus_matrix, *_val ) = phoenix::bind( &PssmAndConsensus::consensus, _1 )
                                ]
                        ]
                    ) > -this->section_separator // sometimes Biobase forget the separator so make it optional
                )
            >   eps[ phoenix::bind( &finalise_matrix, *_val ) ] // tidy up number of sites
            ;
    }
};

/// Partial specialisation to select correct parser
template< typename Iterator >
struct get_entry_parser< Matrix, Iterator > {
    typedef matrix_parser< Iterator > type;
};




/**
 * Parses a fragment entry in a table.
 */
template< typename Iterator >
struct fragment_parser
: qi::grammar< Iterator, Fragment::ptr_t() >
, entry_parser< Iterator >
{
    typedef Fragment parsed_type;
    typedef typename parsed_type::ptr_t parsed_ptr;

    // rules
    qi::rule< Iterator, parsed_ptr() > start;

    fragment_parser() : fragment_parser::base_type( start ) {
        using namespace phoenix;
        using namespace qi;
        using namespace keys;
        using qi::_1;

        // The keys that are used are commented out.
        this->ignored_key.add
//            ( "AC  " ) // Accession number
            ( "AS  " ) // Accession numbers, secondary
            ( "ID  " ) // Identifier
            ( "DT  " ) // Created/Updated
            ( "CO  " ) // Copyright
//            ( "DE  " ) // Description
            ( "OS  " ) // Species
            ( "OC  " ) // Taxonomic classification
            ( "SQ  " ) // Sequence
            ( "SC  " ) // Sequence Source
//            ( "BF  " ) // Binding factors
            ( "PS  " ) // Best supported binding site in the fragment's sequence predicted by Match
            ( "MM  " ) // Method(s)
            ( "DR  " ) // External database links
            ( "RN  " ) // Reference number
            ( "RX  " ) // PUBMED; link to PubMed entry.
            ( "RA  " ) // Reference authors
            ( "RT  " ) // Reference title
            ( "RL  " ) // Reference source
            ;

        start = // parse each section
                eps[ _val = construct< parsed_ptr >( new_< parsed_type >() ) ] // construct the shared_ptr
            >>  +(  (   this->ignored_section // ignore
                    |   this->accession_section[ phoenix::bind( &parsed_type::accession_number, *_val ) = _1 ]
                    |   this->factor_link_section[ phoenix::bind( &parsed_type::factor_links, *_val ) = _1 ]
                    // do not save sequences as data is just too large.
                    // |   this->concatenated_sequence_section [ phoenix::bind( &parsed_type::sequence, *_val ) = _1 ]
                    |   // gene section
                        (   lit( "DE  " )
                        >   (   lit( "Gene: " )
                            >   (   this->tablelink % lit( ", " ) )
                            >   (   lit( " (forward)" )
                                |   lit( " (reverse)" )
                                )
                            ) % lit( "; " )
                        >   eol
                        )
                    ) > -this->section_separator // sometimes Biobase forget the separator so make it optional
                )
            ;
    }
};

/// Partial specialisation to select correct parser
template< typename Iterator >
struct get_entry_parser< Fragment, Iterator > {
    typedef fragment_parser< Iterator > type;
};



/**
 * Parses a site entry in a table.
 */
template< typename Iterator >
struct site_parser
: qi::grammar< Iterator, Site::ptr_t() >
, entry_parser< Iterator >
{
    typedef Site parsed_type;
    typedef typename parsed_type::ptr_t parsed_ptr;

    // rules
    qi::rule< Iterator, parsed_ptr() > start;

    site_parser() : site_parser::base_type( start ) {
        using namespace phoenix;
        using namespace qi;
        using namespace keys;
        using qi::_1;

        // keys that are used are commented out
        this->ignored_key.add
//            ( "AC  " ) // Accession number
            ( "AS  " ) // Accession numbers, secondary
//            ( "ID  " ) // Identifier
            ( "DT  " ) // Created/Updated
            ( "CO  " ) // Copyright
            ( "TY  " ) // Sequence type
//            ( "DE  " ) // Description
            ( "OS  " ) // Species
            ( "OC  " ) // Taxonomic classification
            ( "RE  " ) // Gene region
//            ( "SQ  " ) // Sequence
            ( "EL  " ) // Element
//            ( "S1  " ) // Reference point
//            ( "SF  " ) // Start position
//            ( "ST  " ) // End position
//            ( "BF  " ) // Binding factors
            ( "MX  " ) // Matrices
            ( "SO  " ) // Cellular factor source
            ( "MM  " ) // Method(s)
            ( "CC  " ) // Comments
//            ( "DR  " ) // External database links
            ( "RN  " ) // Reference number
            ( "RX  " ) // PUBMED; link to PubMed entry.
            ( "RA  " ) // Reference authors
            ( "RT  " ) // Reference title
            ( "RL  " ) // Reference source
            ;

        start = // all the sections
                eps[ _val = construct< parsed_ptr >( new_< parsed_type >() ) ]
            >>  +(  (   this->ignored_section // ignore
                    |   this->accession_section[ phoenix::bind( &parsed_type::accession_number, *_val ) = _1 ]
                    |   this->dbref_section[ phoenix::bind( &parsed_type::database_refs, *_val ) = _1 ]
                    |   this->identifier_section[ phoenix::bind( &parsed_type::id, *_val ) = _1 ]
                    |   this->description_section[ phoenix::bind( &parsed_type::description, *_val ) = _1 ]
                    |   this->factor_link_section[ phoenix::bind( &parsed_type::factor_links, *_val ) = _1 ]
                    |   this->multiline_sequence[ phoenix::bind( &parsed_type::sequence, *_val ) = _1 ]
                    |   this->string_line( val( "S1  " ) )[ phoenix::bind( &parsed_type::reference_point, *_val ) = _1 ]
                    |   ( lit( "SF  " ) > int_[ phoenix::bind( &parsed_type::start_position, *_val ) = _1 ] > this->ignored_rest_of_line )
                    |   ( lit( "ST  " ) > int_[ phoenix::bind( &parsed_type::end_position  , *_val ) = _1 ] > this->ignored_rest_of_line )
                    ) > -this->section_separator // sometimes Biobase forget the separator so make it optional
                )
                ;
    }
};

/// Partial specialisation to select correct parser
template< typename Iterator >
struct get_entry_parser< Site, Iterator > {
    typedef site_parser< Iterator > type;
};



/**
 * Parses a biobase table.
 */
template< typename Iterator, typename EntryParser >
struct biobasetable_parser
: qi::grammar< Iterator, typename EntryParser::parsed_type::map_t() >
{
    typedef typename EntryParser::parsed_type parsed_type;
    typedef typename parsed_type::ptr_t parsed_ptr;
    typedef typename EntryParser::parsed_type::map_t map_type;
    typedef typename std::pair< TableLink, parsed_ptr > map_entry;

    // grammars
    EntryParser entry;

    // rules
    qi::rule< Iterator, map_type() > start;
    qi::rule< Iterator, map_entry() > entry_rule;
    qi::rule< Iterator > version_section;

    biobasetable_parser() : biobasetable_parser::base_type( start ) {
        using namespace qi;
        using namespace phoenix;
        using qi::_1;

        version_section =
                +(  lit( "VV  " )
                >   +( char_ - eol )
                >   eol
                )
            >   lit( "XX" )
            >   eol
            >   lit( "//" )
            >   eol
            ;

        entry_rule =
                // make each entry into the value type of the map we will insert it into
                entry[
                    phoenix::bind( &map_entry::first, _val )
                        = phoenix::bind( &parsed_type::accession_number, *_1 ),
                    phoenix::bind( &map_entry::second, _val ) = _1
                ]
            >   lit( "//" )
            >   eol
            ;

        start =
                version_section
            >   +entry_rule
            ;
    }
};





} //namespace spirit
BIO_NS_END




#endif //BIO_JR_30AUG2011_TRANSFAC_QI_H_

