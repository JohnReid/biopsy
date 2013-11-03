/** Copyright John Reid 2012
 *
 * \file
 * \brief Outputs TRANSFAC data structures using boost::spirit karma.
 */


#ifndef BIO_JR_6SEP2011_TRANSFAC_KARMA_H_
#define BIO_JR_6SEP2011_TRANSFAC_KARMA_H_



#include "bio/spirit/transfac_fusion.h"

#include <boost/spirit/include/karma.hpp>
#include <boost/spirit/include/support_adapt_adt_attributes.hpp>




namespace boost {
namespace spirit {
namespace traits {

// specialise how iterators into containers of pointers are dereferenced
template <>
struct deref_iterator< karma::detail::indirect_iterator< BIO_NS::AlignDescList::const_iterator > >
{
    typedef BIO_NS::AlignDesc const & type;
    typedef karma::detail::indirect_iterator< BIO_NS::AlignDescList::const_iterator > iterator;

    static
    type
    call( iterator const & it ) {
    	return **it;
    }
};


// specialise how iterators into containers of pointers are dereferenced
template <>
struct deref_iterator< karma::detail::indirect_iterator< BIO_NS::FactorLinkList::const_iterator > >
{
    typedef BIO_NS::FactorLink const & type;
    typedef karma::detail::indirect_iterator< BIO_NS::FactorLinkList::const_iterator > iterator;

    static
    type
    call( iterator const & it ) {
    	return **it;
    }
};


} // namespace traits
} // namespace spirit
} // namespace boost


BIO_NS_START
namespace spirit {

namespace karma = boost::spirit::karma;
namespace ascii = boost::spirit::ascii;
namespace phoenix = boost::phoenix;



/// Specialised to select correct generator for Entry type.
template< typename Entry, typename Iterator >
struct get_entry_generator;



/**
 * Generates a binding factor.
 */
template< typename Iterator >
struct factorlink_generator
: karma::grammar< Iterator, FactorLink() >
{
    karma::rule< Iterator, FactorLink() > start;

    factorlink_generator() : factorlink_generator::base_type( start )
    {
        using namespace phoenix;
        using namespace karma;
        using karma::_1;

        start =
        	stream // the table link
        	<< lit("; ") << *char_ // the name
        	<< lit("; Quality: ") << int_
        	<< lit("; Cellular source: ") << *char_
        	<< lit("; Sites included: ") << bool_
        	<< lit("; Species: ") << +char_ % lit(',')
        	<< lit(".")
        	;
    }
};



/**
 * Generates a binding site.
 */
template< typename Iterator >
struct aligndesc_generator
: karma::grammar< Iterator, AlignDesc() >
{
    karma::rule< Iterator, AlignDesc() > start;
	karma::rule< Iterator, std::vector<int> const &() > gaps;

    aligndesc_generator() : aligndesc_generator::base_type( start )
    {
        using namespace karma;
        using karma::_1;

        gaps = int_ % lit(',');
        start =
        	*char_ // the sequence
        	<< lit("; ") << stream // the link
        	<< skip[lit(", ") << stream] // the secondary link
        	<< lit("; ") << int_ // the start
        	<< lit("; ") << int_ // the length
        	<< lit("; ") << -gaps // the gaps
        	<< lit("; ") << bool_ // positive orientation
        	<< lit(".")
			;
    }
};



/**
 * Generates a PSSM entry.
 */
template< typename Iterator >
struct pssmentry_generator
: karma::grammar< Iterator, PssmEntry() >
{
    karma::rule< Iterator, PssmEntry() > start;

    pssmentry_generator() : pssmentry_generator::base_type( start )
    {
        using namespace karma;

        start =
        	float_ << lit("  ")
        	<< float_ << lit("  ")
        	<< float_ << lit("  ")
        	<< float_
        	;
    }
};



/**
 * Base class for generators of biobase table entries.
 */
template< typename Iterator >
struct entry_generator
{
	// rules
    karma::rule< Iterator, std::string() > string;
    karma::rule< Iterator, TableLink() > tablelink;
    karma::rule< Iterator, TableLinkVec() > tablelinkvec;
    karma::rule< Iterator, Identifier() > identifier;
    karma::rule< Iterator, FactorLink() > factorlink_line;
    karma::rule< Iterator, AlignDesc() > aligndesc_line;
    karma::rule< Iterator, db_ref() > dbref_line;

    // generators
    factorlink_generator< Iterator > factorlink;
    aligndesc_generator< Iterator > aligndesc;

    entry_generator()
    {
    	using namespace karma;

        string             = *char_;
        tablelink          = stream;
        tablelinkvec       = -( this->tablelink % "; " );
        identifier         = stream;
        factorlink_line    = lit("BF  ") << factorlink << eol;
        aligndesc_line     = lit("BS  ") << aligndesc << eol;
        dbref_line         = lit("DR  ") << stream << eol;
    }
};



/**
 * Generates a matrix.
 */
template< typename Iterator >
struct matrix_generator
: karma::grammar< Iterator, Matrix() >
, entry_generator< Iterator >
{
	typedef Matrix generated_type;

	// rules
    karma::rule< Iterator, generated_type() > start;
    karma::rule< Iterator, Pssm(), karma::locals< int > > pssm;
    karma::rule< Iterator, ConsensusMatrix const &() > consensusmatrix;

    // generators
    pssmentry_generator< Iterator > pssmentry;

    matrix_generator() : matrix_generator::base_type( start )
    {
        using namespace karma;

        pssm =
        	//eps[ _a = 0 ]
        	lit("P0      A      C      G      T") << eol
        	<< *( lit("--  ") << pssmentry << eol )
        	//<< *pssmentry
        	;
        consensusmatrix = *char_;
        start =
        	lit("AC  ") << this->tablelink << eol
        	<< lit("ID  ") << this->identifier << eol
        	<< lit("NA  ") << this->string << eol
        	<< lit("BA  ") << this->string << eol
        	<< lit("--  Number of sites: ") << int_ << eol
        	<< lit("DE  ") << this->string << eol
        	<< *this->factorlink_line
        	<< *this->aligndesc_line
        	<< pssm
        	<< lit("--  Consensus: ") << consensusmatrix << eol
        	;
    }
};


/// Partial specialisation to select correct generator
template< typename Iterator >
struct get_entry_generator< Matrix, Iterator > {
	typedef matrix_generator< Iterator > type;
};




/**
 * Generates a site.
 */
template< typename Iterator >
struct site_generator
: karma::grammar< Iterator, Site() >
, entry_generator< Iterator >
{
	typedef Site generated_type;

	// rules
    karma::rule< Iterator, generated_type() > start;

    // generators

    site_generator() : site_generator::base_type( start )
    {
        using namespace karma;

        start =
        	lit("AC  ") << this->tablelink << eol
        	<< lit("ID  ") << this->identifier << eol
        	<< lit("SQ  ") << this->string << eol
        	<< lit("DE  ") << this->string << eol
        	<< *this->factorlink_line
        	<< *this->dbref_line
        	<< lit("S1  ") << this->string << eol
        	<< lit("SF  ") << int_ << eol
        	<< lit("ST  ") << int_ << eol
        	;
    }
};


/// Partial specialisation to select correct generator
template< typename Iterator >
struct get_entry_generator< Site, Iterator > {
	typedef site_generator< Iterator > type;
};




/**
 * Generates a factor.
 */
template< typename Iterator >
struct factor_generator
: karma::grammar< Iterator, Factor() >
, entry_generator< Iterator >
{
    typedef Factor generated_type;

    karma::rule< Iterator, generated_type() > start;

    factor_generator() : factor_generator::base_type( start )
    {
        using namespace karma;

        start =
                lit("AC  ") << this->tablelink << eol
            <<  lit("FA  ") << this->string << eol
            <<  lit("TY  ") << stream << eol
            <<  lit("GE  ") << this->tablelink << eol
            <<  lit("MX  ") << this->tablelinkvec << eol
            <<  lit("SY  ") << -( this->string % "; " ) << eol
            <<  *this->dbref_line
            <<  lit("OC  ") << -( this->string % "; " ) << eol
            <<  lit("ST  ") << this->tablelinkvec << eol
            <<  lit("CX  ") << this->tablelinkvec << eol
            <<  lit("HC  ") << this->tablelinkvec << eol
            <<  lit("HP  ") << this->tablelinkvec << eol
            ;
    }
};


/// Partial specialisation to select correct generator
template< typename Iterator >
struct get_entry_generator< Factor, Iterator > {
    typedef factor_generator< Iterator > type;
};


/**
 * Generates a fragment.
 */
template< typename Iterator >
struct fragment_generator
: karma::grammar< Iterator, Fragment() >
, entry_generator< Iterator >
{
    typedef Fragment generated_type;

    karma::rule< Iterator, generated_type() > start;

    fragment_generator() : fragment_generator::base_type( start )
    {
        using namespace karma;

        start =
                lit("AC  ") << this->tablelink << eol
            <<  *this->factorlink_line
            <<  lit("GE  ") << this->tablelinkvec << eol
            <<  lit("SQ  ") << this->string << eol
            ;
    }
};


/// Partial specialisation to select correct generator
template< typename Iterator >
struct get_entry_generator< Fragment, Iterator > {
    typedef fragment_generator< Iterator > type;
};


/**
 * Generates a gene.
 */
template< typename Iterator >
struct gene_generator
: karma::grammar< Iterator, Gene() >
, entry_generator< Iterator >
{
    typedef Gene generated_type;

    karma::rule< Iterator, generated_type() > start;

    gene_generator() : gene_generator::base_type( start )
    {
        using namespace karma;

        start =
                lit("AC  ") << this->tablelink << eol
            <<  lit("ID  ") << this->string << '$' << this->string << eol
            <<  *this->dbref_line
            ;
    }
};


/// Partial specialisation to select correct generator
template< typename Iterator >
struct get_entry_generator< Gene, Iterator > {
    typedef gene_generator< Iterator > type;
};



/**
 * Generates an entry from the compel table.
 */
template< typename Iterator >
struct compel_generator
: karma::grammar< Iterator, Compel() >
, entry_generator< Iterator >
{
    typedef Compel generated_type;

    karma::rule< Iterator, generated_type() > start;

    compel_generator() : compel_generator::base_type( start )
    {
        using namespace karma;

        start =
                lit("AC  ") << this->tablelink << eol
            <<  lit("ID  ") << this->identifier << eol
            <<  lit("GE  ") << this->tablelink << eol
            <<  lit("SQ  ") << this->string << eol
            <<  lit("PS  ") << int_ << " to " << int_ << eol
            //<<  *this->dbref_line
            ;
    }
};


/// Partial specialisation to select correct generator
template< typename Iterator >
struct get_entry_generator< Compel, Iterator > {
    typedef compel_generator< Iterator > type;
};


/**
 * Generates an entry from the evidence table.
 */
template< typename Iterator >
struct evidence_generator
: karma::grammar< Iterator, Evidence() >
, entry_generator< Iterator >
{
    typedef Evidence generated_type;

    karma::rule< Iterator, generated_type() > start;

    evidence_generator() : evidence_generator::base_type( start )
    {
        using namespace karma;

        start =
                lit("AC  ") << this->tablelink << eol
            <<  lit("CE  ") << this->tablelink << eol
            <<  *this->dbref_line
            ;
    }
};


/// Partial specialisation to select correct generator
template< typename Iterator >
struct get_entry_generator< Evidence, Iterator > {
    typedef evidence_generator< Iterator > type;
};



/**
 * Generates an entry from the Molecule table.
 */
template< typename Iterator >
struct molecule_generator
: karma::grammar< Iterator, Molecule() >
, entry_generator< Iterator >
{
    typedef Molecule generated_type;

    karma::rule< Iterator, generated_type() > start;

    molecule_generator() : molecule_generator::base_type( start )
    {
        using namespace karma;

//        start =
//                lit("AC  ") << this->tablelink << eol
//            ;
    }
};


/// Partial specialisation to select correct generator
template< typename Iterator >
struct get_entry_generator< Molecule, Iterator > {
    typedef molecule_generator< Iterator > type;
};



/**
 * Generates an entry from the Pathway table.
 */
template< typename Iterator >
struct pathway_generator
: karma::grammar< Iterator, Pathway() >
, entry_generator< Iterator >
{
    typedef Pathway generated_type;

    karma::rule< Iterator, generated_type() > start;

    pathway_generator() : pathway_generator::base_type( start )
    {
        using namespace karma;

//        start =
//                lit("AC  ") << this->tablelink << eol
//            ;
    }
};


/// Partial specialisation to select correct generator
template< typename Iterator >
struct get_entry_generator< Pathway, Iterator > {
    typedef pathway_generator< Iterator > type;
};



} //namespace spirit
BIO_NS_END

#endif //BIO_JR_6SEP2011_TRANSFAC_KARMA_H_

