/**
@file

Copyright John Reid 2006

*/
#ifdef _MSC_VER
# pragma warning( disable : 4503 )
#endif

#include "biopsy/custom_pssm.h"
#include "biopsy/sequence.h"
#include "biopsy/transfac.h"

#include <bio/biobase_filter.h>
#include <bio/cache.h>
#include <bio/singleton.h>
#include <bio/environment.h>

#include <boost/filesystem/convenience.hpp>
#include <boost/spirit/home/classic/phoenix/binders.hpp>
#include <boost/spirit/home/classic/core/composite/kleene_star.hpp>
#include <boost/spirit/home/classic/utility/lists.hpp>
#include <boost/spirit/include/classic_file_iterator.hpp>
#include <boost/spirit/include/classic_core.hpp>
#include <boost/spirit/include/classic_attribute.hpp>
#if BOOST_VERSION >= 103500
 using namespace BOOST_SPIRIT_CLASSIC_NS;
#endif

using namespace boost::spirit;
using namespace phoenix;
using namespace boost;
using namespace std;


#ifdef WIN32
# define DIR_SEP "\\"
#else //WIN32
# define DIR_SEP "/"
#endif //WIN32



namespace biopsy {

namespace detail {

    
struct CustomPssmFromName
    : std::unary_function< std::string, custom_pssm::ptr >
{
    custom_pssm::ptr operator()( const std::string & name ) const
    {
        return parse_custom_pssm_file( custom_pssm_filename( name ) ); //the name should be the PSSM id
    }
};

struct custom_pssm_cache
    : BIO_NS::Cache< CustomPssmFromName >
    , BIO_NS::Singleton< custom_pssm_cache >
{
};


struct custom_pssm_info
{
    string name;
    string id;
    string url;
    vector< unsigned > positions;
    vector< double > counts_a, counts_c, counts_g, counts_t;
    vector< string > iupac;
    unsigned width;

    custom_pssm::ptr create_pssm() const
    {
        //calculate pssm
        custom_pssm::ptr result( new custom_pssm );
        result->name = name;
        result->id = id;
        result->url = url;
        result->counts.resize( width );
        for( unsigned i = 0; width != i; ++i )
        {
            //const double sum = counts_a[i] + counts_c[i] + counts_g[i] + counts_t[i];
            result->counts[i].set( 0, counts_a[i] );
            result->counts[i].set( 1, counts_c[i] );
            result->counts[i].set( 2, counts_g[i] );
            result->counts[i].set( 3, counts_t[i] );
        }
        return result;
    }
};

struct custom_pssm_closure : BOOST_SPIRIT_CLASSIC_NS::closure< custom_pssm_closure, custom_pssm_info >
{
    member1 pssm;
};

struct push_back_impl {

    template< typename Arg1, typename Arg2 >
    struct result { typedef void type; };

    template< typename Arg1, typename Arg2 >
    void operator()( Arg1 & container, Arg2 value ) const
    { 
        container.push_back(value);
    }
};

struct insert_impl {

    template< typename Arg1, typename Arg2 >
    struct result { typedef void type; };

    template< typename Arg1, typename Arg2 >
    void operator()( Arg1 & container, Arg2 value ) const
    { 
        container.insert(value);
    }
};

struct custom_pssm_grammar : public grammar< custom_pssm_grammar, custom_pssm_closure::context_t >
{
    template< typename ScannerT >
    struct definition
    {
        typedef rule< ScannerT > rule_t;
        rule_t custom_pssm_rule, line_rule, counts_a_rule, counts_c_rule, counts_g_rule, counts_t_rule, position_rule, iupac_rule;
        definition( custom_pssm_grammar const& g )
        { 
            using phoenix::function;

            member_var_ptr< string, custom_pssm_info > name = &custom_pssm_info::name;
            member_var_ptr< string, custom_pssm_info > id = &custom_pssm_info::id;
            member_var_ptr< string, custom_pssm_info > url = &custom_pssm_info::url;
            member_var_ptr< unsigned, custom_pssm_info > width = &custom_pssm_info::width;
            member_var_ptr< vector< double >, custom_pssm_info > counts_a = &custom_pssm_info::counts_a;
            member_var_ptr< vector< double >, custom_pssm_info > counts_c = &custom_pssm_info::counts_c;
            member_var_ptr< vector< double >, custom_pssm_info > counts_g = &custom_pssm_info::counts_g;
            member_var_ptr< vector< double >, custom_pssm_info > counts_t = &custom_pssm_info::counts_t;
            member_var_ptr< vector< unsigned >, custom_pssm_info > positions = &custom_pssm_info::positions;
            member_var_ptr< vector< string >, custom_pssm_info > iupac = &custom_pssm_info::iupac;

            function< push_back_impl > push_back;

            counts_a_rule = list_p.direct( real_p[ push_back( counts_a(g.pssm), arg1 ) ], eps_p );
            counts_c_rule = list_p.direct( real_p[ push_back( counts_c(g.pssm), arg1 ) ], eps_p );
            counts_g_rule = list_p.direct( real_p[ push_back( counts_g(g.pssm), arg1 ) ], eps_p );
            counts_t_rule = list_p.direct( real_p[ push_back( counts_t(g.pssm), arg1 ) ], eps_p );
            position_rule = list_p.direct( uint_p[ push_back( positions(g.pssm), arg1 ) ], eps_p );
            iupac_rule = list_p.direct( (+range_p('A','Z'))[ push_back( iupac(g.pssm), construct_< string >(arg1, arg2) ) ], eps_p );

            line_rule =
                ( ch_p('#') >> *(print_p-eol_p) )
                | ( str_p("NA  ") >> (+(print_p-eol_p))[ assign_a( name(g.pssm)() ) ] )
                | ( str_p("ID  ") >> (+(print_p-eol_p))[ assign_a( id(g.pssm)() ) ] )
                | ( str_p("WI  ") >> uint_p[ width(g.pssm) = arg1 ] )
                | ( str_p("CA  ") >> counts_a_rule )
                | ( str_p("CC  ") >> counts_c_rule )
                | ( str_p("CG  ") >> counts_g_rule )
                | ( str_p("CT  ") >> counts_t_rule )
                | ( str_p("PO  ") >> position_rule )
                | ( str_p("IU  ") >> iupac_rule )
                | ( str_p("UR  ") >> (+graph_p)[ assign_a( url(g.pssm)() ) ] )
                | eps_p
                ;

            custom_pssm_rule = list_p( line_rule, +eol_p ) ;
        }
        rule_t const& start() const { return custom_pssm_rule; }
    };
};


struct pssm_set_closure : BOOST_SPIRIT_CLASSIC_NS::closure< pssm_set_closure, pssm_set >
{
    member1 _set;
};

struct pssm_set_grammar : public grammar< pssm_set_grammar, pssm_set_closure::context_t >
{
    template< typename ScannerT >
    struct definition
    {
        typedef rule< ScannerT > rule_t;
        rule_t pssm_set_rule, line_rule;
        definition( pssm_set_grammar const& g )
        { 
            using phoenix::function;

            function< insert_impl > insert;

            line_rule =
                ( ch_p('#') >> *(print_p-eol_p) ) //ignore comments
                | list_p.direct( (+graph_p)[ insert( g._set, construct_< std::string >( arg1, arg2 ) ) ], blank_p ) //whitespace separated list of strings
                | eps_p //empty line
                ;

            pssm_set_rule = list_p( line_rule, +eol_p ) ;
        }
        rule_t const& start() const { return pssm_set_rule; }
    };
};

} //namespace detail



custom_pssm::ptr parse_custom_pssm_file( const std::string & custom_pssm_filename )
{
    using namespace boost::spirit;
    using namespace std;

    //open the file
    file_iterator<> first( custom_pssm_filename.c_str() );
    if( !first ) throw std::logic_error(BIOPSY_MAKE_STRING( "Could not open: " << custom_pssm_filename ));
    file_iterator<> last = first.make_end();

    typedef char                     char_t;
    typedef file_iterator<char_t>    iterator_t;

    //parse the file
    detail::custom_pssm_grammar g;
    detail::custom_pssm_info pssm_info;
    parse_info<iterator_t> info = parse( first, last, g[ var(pssm_info)=arg1 ], blank_p );
    if( ! info.full ) 
    {
        const string stopped_at(info.stop, last);
        throw std::logic_error(BIOPSY_MAKE_STRING( "Could not parse entire file, stopped at: \"" << stopped_at << "\""));    
    }

    //done parsing - check what we have makes sense
    if( pssm_info.width != pssm_info.counts_a.size()
        || pssm_info.width != pssm_info.counts_c.size()
        || pssm_info.width != pssm_info.counts_g.size()
        || pssm_info.width != pssm_info.counts_t.size()
        //|| pssm_info.counts_a.size() != pssm_info.positions.size()
        //|| pssm_info.counts_a.size() != pssm_info.iupac.size()
    ) throw std::logic_error( "Wrong number of entries in A, C, G, or T counts" );

    return pssm_info.create_pssm();
}



std::string custom_pssm_filename( const std::string & pssm_id )
{
    return BIOPSY_MAKE_STRING( 
        BIO_NS::BioEnvironment::singleton().get_custom_pssm_dir()
        << DIR_SEP
        << pssm_id
        << ".pssm" )
        ;
}


custom_pssm::ptr get_custom_pssm( const std::string & pssm_id )
{
    return detail::custom_pssm_cache::singleton()( pssm_id );
}

pssm_set_ptr get_pssm_set( const std::string & name )
{
    namespace fs = boost::filesystem;
    const fs::path custom_pssm_dir = fs::system_complete( fs::path( BIO_NS::BioEnvironment::singleton().get_custom_pssm_dir() ) );
    const fs::path pssm_set_file = custom_pssm_dir / (name + ".pssm_set");

    //open the file
    file_iterator<> first( pssm_set_file.string() );
    if( !first ) throw std::logic_error(BIOPSY_MAKE_STRING( "Could not open: " << pssm_set_file.string() ));
    file_iterator<> last = first.make_end();

    typedef char                     char_t;
    typedef file_iterator<char_t>    iterator_t;

    //parse the file
    detail::pssm_set_grammar g;
    pssm_set _pssm_set;
    parse_info< iterator_t > info = parse( first, last, g[ var(_pssm_set) = arg1 ] );
    if( ! info.full ) 
    {
        const string stopped_at(info.stop, last);
        throw std::logic_error(BIOPSY_MAKE_STRING( "Could not parse entire file, stopped at: \"" << stopped_at << "\""));    
    }

    pssm_set_ptr result( new pssm_set( _pssm_set ) );

    return result;
}

pssm_set_ptr all_custom_pssms()
{
    USING_BIO_NS;

    pssm_set_ptr result( new pssm_set );

    namespace fs = boost::filesystem;
    const std::string custom_pssm_dir = BioEnvironment::singleton().get_custom_pssm_dir();
    const fs::path full_path = fs::system_complete( fs::path( custom_pssm_dir ) );
    if ( !fs::exists( full_path ) ) throw std::logic_error( BIOPSY_MAKE_STRING( "Custom pssm dir does not exist: " << full_path.string() ) );
    if ( ! fs::is_directory( full_path ) ) throw std::logic_error( BIOPSY_MAKE_STRING( "Custom pssm dir is not a directory: " << full_path.string() ) );
    fs::directory_iterator end_iter;
    for( fs::directory_iterator dir_itr( full_path );
          dir_itr != end_iter;
          ++dir_itr )
    {
        if ( fs::extension( dir_itr->path() ) == ".pssm" )
        {
            result->insert( fs::basename( dir_itr->path() ) );
        }
    }

    return result;
}


std::vector< std::string > get_installed_custom_pssm_set_names()
{
    //implementation essentially same as all_custom_pssms()
    USING_BIO_NS;

    std::vector< std::string > result;

    namespace fs = boost::filesystem;
    const std::string custom_pssm_dir = BioEnvironment::singleton().get_custom_pssm_dir();
    const fs::path full_path = fs::system_complete( fs::path( custom_pssm_dir ) );
    if ( !fs::exists( full_path ) ) throw std::logic_error( BIOPSY_MAKE_STRING( "Custom pssm dir does not exist: " << full_path.string() ) );
    if ( ! fs::is_directory( full_path ) ) throw std::logic_error( BIOPSY_MAKE_STRING( "Custom pssm dir is not a directory: " << full_path.string() ) );
    fs::directory_iterator end_iter;
    for( fs::directory_iterator dir_itr( full_path );
          dir_itr != end_iter;
          ++dir_itr )
    {
        if ( fs::extension( dir_itr->path() ) == ".pssm_set" )
        {
            result.push_back( fs::basename( dir_itr->path() ) );
        }
    }

    return result;
}


void add_pssms_from_pssm_set(
    pssm_set & pssms,
    const std::string & pssm_set_name,
    bool use_transfac_consensus_sequences,
    const std::string & matrix_species,
    const std::string & matrix_name_match )
{
    if( pssm_set_name == "transfac" )
    {
        BIO_NS::BiobasePssmFilter filter(
            use_transfac_consensus_sequences,
            matrix_species,
            matrix_name_match );
        string_vec_ptr transfac_pssm_names = get_transfac_pssm_accessions( filter );
        std::copy( transfac_pssm_names->begin(), transfac_pssm_names->end(), inserter( pssms, pssms.end() ) );
    }
    else
    {
        pssm_set_ptr _pssm_set = pssm_set_name == "all-custom" ? all_custom_pssms() : get_pssm_set( pssm_set_name );
        std::copy( _pssm_set->begin(), _pssm_set->end(), inserter( pssms, pssms.end() ) );
    }
}


/**
Get the list of the custom pssms using boost python compatible type*/
string_vec_ptr
get_all_custom_pssms()
{
    pssm_set_ptr all_pssms = all_custom_pssms();
    string_vec_ptr result( new string_vec );
    BOOST_FOREACH(std::string pssm, *all_pssms)
        result->push_back(pssm);
    return result;
}

/**
Get the list of the custom pssms in boost python compatible type*/
string_vec_ptr
get_custom_pssms(const std::string & pssm_set_name)
{
    pssm_set_ptr all_pssms = get_pssm_set(pssm_set_name);
    string_vec_ptr result( new string_vec );
    BOOST_FOREACH(std::string pssm, *all_pssms)
        result->push_back(pssm);
    return result;
}

/**
Get the list of the custom pssms in boost python compatible type*/
string_vec_ptr
get_custom_pssm_set_names()
{
    std::vector< std::string > names = get_installed_custom_pssm_set_names();

    string_vec_ptr result( new string_vec );
    BOOST_FOREACH(std::string name,names)
    {
        result->push_back(name);
    }
    return result;
}



} //namespace biopsy
