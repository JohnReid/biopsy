/**
@file

Copyright John Reid 2006-2011

*/

#include "biopsy/pssm.h"
#include "biopsy/custom_pssm.h"
#include "biopsy/sequence.h"

#include <bio/biobase_db.h>
#include <bio/sequence.h>
#include <bio/cache.h>
#include <bio/singleton.h>
#include <bio/serialisable.h>
#include <bio/environment.h>

using namespace boost;
using namespace std;



namespace biopsy {


pssm_parameters::pssm_parameters()
    : pseudo_counts( 0.25 )
    , likelihoods_size( 100 )
    , calculate_likelihoods_map_size( 10000 )
    , binding_background_odds_prior( 2e-5 )
    , use_cumulative_dists( true )
    , use_p_value( false )
    , use_score( false )
    , avg_phylo_bayes( true )
    , max_chain_num_boxes_limit( 50000 )
    , min_related_evidence_fraction( .5 )
{
}

nucleo_dist::nucleo_dist(
    double a,
    double c,
    double g,
    double t )
{
    _values[ 0 ] = a;
    _values[ 1 ] = c;
    _values[ 2 ] = g;
    _values[ 3 ] = t;
}

nucleo_dist::vec  operator + (const nucleo_dist::vec & v,double pseudo_count)
{
    nucleo_dist::vec retVal;

    BOOST_FOREACH(const nucleo_dist & c, v)
    {
        nucleo_dist dist(
            c.get(0) + pseudo_count,
            c.get(1) + pseudo_count,
            c.get(2) + pseudo_count,
            c.get(3) + pseudo_count);
        retVal.push_back(dist);
    }
    return retVal;
}



/** argument is 0 for a, 1 for c, ... */
double
nucleo_dist::get( unsigned nucleo ) const
{
    return _values[ nucleo ];
}

void
nucleo_dist::set( unsigned nucleo, double value )
{
    _values[ nucleo ] = value;
}

double
nucleo_dist::get_total( ) const
{
    return std::accumulate( _values.begin(), _values.end(), 0.0 );
}

/** argument is 0 for a, 1 for c, ... */
double
nucleo_dist::get_freq( unsigned nucleo ) const
{
    return get( nucleo ) / get_total();
}

bool
nucleo_dist::operator==( const nucleo_dist & rhs ) const
{
    return
        _values == rhs._values;
}

nucleo_dist
uniform_nucleo_dist()
{
    return nucleo_dist( 0.25, 0.25, 0.25, 0.25 );
}

nucleo_dist
get_dist_for_iupac( char s )
{
    USING_BIO_NS;
    IupacCode i( s );
    return
        nucleo_dist(
            i.get_num( 'a' ),
            i.get_num( 'c' ),
            i.get_num( 'g' ),
            i.get_num( 't' ) );
}



pssm_ptr
create_pssm(
    const nucleo_dist::vec & dists )
{
    //from http://nar.oxfordjournals.org/cgi/reprint/31/13/3576

    //first calculate the information vector and the min and max scores possible
    std::vector< double > infos;
    std::vector< double > min_freqs;
    double min_score = 0.0;
    double max_score = 0.0;
    BOOST_FOREACH( const nucleo_dist & dist, dists )
    {
        double info = 0.0;
        double min_freq = 1.0;
        double max_freq = 0.0;
        for( unsigned i = 0; 4 != i; ++i )
        {
            const double fib = dist.get_freq( i );
            if( 0.0 > fib || fib > 1.0 ) {
                throw std::logic_error( BIOPSY_MAKE_STRING( "pssm_create(): Frequency is not in [0,1]: "<<fib ) );
            }
            min_freq = std::min( fib, min_freq );
            max_freq = std::max( fib, max_freq );
            if( 0.0 != fib )
            {
                info += fib * std::log( 4.0 * fib );
            }
        }
        //std::cout << info << "\n";
        min_score += min_freq * info;
        max_score += max_freq * info;
        if( info < 0. ) {
            throw std::logic_error( BIOPSY_MAKE_STRING( "pssm_create(): info < 0.0: "<<info ) );
        }
        infos.push_back( info );
        if( min_freq < 0. ) {
            throw std::logic_error( BIOPSY_MAKE_STRING( "pssm_create(): min_freq < 0.0: "<<info ) );
        }
        min_freqs.push_back( min_freq );
    }

    if( max_score == min_score ) {
        throw std::logic_error( BIOPSY_MAKE_STRING( "pssm_create(): max score == min score: "<<min_score ) );
    }
    if( max_score < min_score ) {
        throw std::logic_error( BIOPSY_MAKE_STRING( "pssm_create(): max score < min score: "<<max_score<<" < "<<min_score ) );
    }
    const double range = 1.0 / ( max_score - min_score );
    if( BIOPSY_ISNAN( range ) ) {
        throw std::logic_error( "pssm_create(): Range is too small to calculate" );
    }

    //now calculate the scores for each position
    pssm_ptr result( new pssm );
    for( unsigned pos = 0; dists.size() != pos; ++pos )
    {
        const nucleo_dist scores(
            infos[ pos ] * ( std::max( 0.0, dists[ pos ].get_freq( 0 ) - min_freqs[ pos ] ) ) * range,
            infos[ pos ] * ( std::max( 0.0, dists[ pos ].get_freq( 1 ) - min_freqs[ pos ] ) ) * range,
            infos[ pos ] * ( std::max( 0.0, dists[ pos ].get_freq( 2 ) - min_freqs[ pos ] ) ) * range,
            infos[ pos ] * ( std::max( 0.0, dists[ pos ].get_freq( 3 ) - min_freqs[ pos ] ) ) * range );
        for( unsigned b = 0; 4 != b; ++b ) {
            if( scores.get( b ) < 0.0 ) {
                throw std::logic_error( BIOPSY_MAKE_STRING( "pssm_create(): Have a negative score in PSSM: "<<scores.get( b ) ) );
            }
        }
        result->push_back( scores );
    }

    return result;
}

inline
unsigned
get_nucleo_index( const sequence::value_type & c )
{
    switch( c )
    {
    case 'a':
    case 'A':
        return 0;
    case 'c':
    case 'C':
        return 1;
    case 'g':
    case 'G':
        return 2;
    case 't':
    case 'T':
        return 3;
    }

    throw std::logic_error( BIOPSY_MAKE_STRING( "Not a nucleotide: \"" << c << "\"" ) );
}


double
score(
    const pssm & pssm,
    sequence::const_iterator s_begin )
{
    double score = 0.0;
    for( unsigned pos = 0; pssm.size() != pos; ++pos, ++s_begin )
    {
        score += ( pssm [ pos ] ).get( get_nucleo_index( *s_begin ) );
    }

    return score;
}

double
score_complement(
    const pssm & pssm,
    sequence::const_iterator s_begin )
{
    double score = 0.0;
    const unsigned size = pssm.size();
    for( unsigned pos = 0; size != pos; ++pos, ++s_begin )
    {
        score += ( pssm [ size - pos - 1 ] ).get( get_nucleo_index( nucleo_complement()( *s_begin ) ) );
    }

    return score;
}

double
score_pssm(
    pssm_ptr pssm,
    const sequence & s )
{
    if( pssm->size() > s.size() )
    {
        throw std::logic_error( "Pssm too large for sequence" );
    }

    return score( *pssm, s.begin() );
}


double
get_odds_ratio(
    double p_score_under_binding,
    double p_score_under_background )
{
    return pssm_parameters::singleton().binding_background_odds_prior * p_score_under_binding / p_score_under_background;
}



double
get_odds_ratio(
    double score,
    likelihoods_ptr _binding_dist,
    likelihoods_ptr _background_dist )
{
    const double p_score_under_binding = get_likelihood( _binding_dist, score );
    const double p_score_under_background = get_likelihood( _background_dist, score );

    return get_odds_ratio( p_score_under_binding, p_score_under_background );
}


double
get_p_binding(
    double odds_ratio )
{
    return odds_ratio / ( 1.0 + odds_ratio );
}


double
get_odds_ratio_from_p_binding(
    double p_binding )
{
    return p_binding / ( 1.0 - p_binding );
}


namespace detail {

enum pssm_type {
    PSSM_TRANSFAC,
    PSSM_CUSTOM,
    PSSM_UNKNOWN_TYPE
};

pssm_type get_pssm_type( const std::string & id );


struct PssmFromName
    : std::unary_function< std::string, pssm_info >
{
    const static regex transfac_re; /**< Regex to match TRANSFAC PSSM names. */
    const static regex custom_re; /**< Regex to match custom PSSM names. */
    pssm_info operator()( const std::string & name ) const
    {
        std::cout << "Creating pssm info for: " << name << "\n";

        nucleo_dist::vec counts;
        int n_sites = 0;

        switch( get_pssm_type( name ) )
        {
        case PSSM_TRANSFAC:
            {
                //look in transfac
                USING_BIO_NS;
                const TableLink link = parse_table_link_accession_number( name );

                switch( link.table_id )
                {
                case SITE_DATA:
                    {
                        Site::map_t::const_iterator s = BiobaseDb::singleton().get_sites().find( link );
                        if( BiobaseDb::singleton().get_sites().end() == s )
                        {
                            throw std::logic_error( BIOPSY_MAKE_STRING( "Could not find TRANSFAC PSSM id: " << name ) );
                        }
                        //std::cout << s->second->sequence << "\n";
                        BOOST_FOREACH( char c, s->second->sequence )
                        {
                            const nucleo_dist count = get_dist_for_iupac( c );
                            counts.push_back( count );
                        }
                        //    These are consensus sequences, so assume they are derived from 
                        //    12 sites
                        n_sites = 12;
                    }
                    break;

                case MATRIX_DATA:
                    {
                        Matrix::map_t::const_iterator m = BiobaseDb::singleton().get_matrices().find( link );
                        if( BiobaseDb::singleton().get_matrices().end() == m )
                        {
                            throw std::logic_error( BIOPSY_MAKE_STRING( "Could not find TRANSFAC PSSM id: " << name ) );
                        }
                        BOOST_FOREACH( const PssmEntry & e, m->second->pssm )
                        {
                            const nucleo_dist count(
                                e.get_count( 'a' ),
                                e.get_count( 'c' ),
                                e.get_count( 'g' ),
                                e.get_count( 't' ) );
                            counts.push_back( count );
                        }
                        n_sites =  m->second->number_of_sites;
                    }
                    break;

                default:
                    throw std::logic_error( BIOPSY_MAKE_STRING( "Could not parse TRANSFAC PSSM id: " << name ) );
                }
            }
            break;

        case PSSM_CUSTOM:
            {
                //it looks like the name of a custom PSSM
                custom_pssm::ptr pssm = parse_custom_pssm_file( custom_pssm_filename( name ) ); //the name should be the PSSM id
                counts = pssm->counts;
                n_sites = 0;
            }
            break;

        default:
            throw std::logic_error( BIOPSY_MAKE_STRING( "Unknown PSSM name: " << name ) );
        }

        //add pseudo-counts
        double max_count = 0;
        double min_count = 999;
        double multiplier = 1;

        //    For distributions where the sum adds up to 1, divide the pseudo count by 100;
        //    Have to check for total less than 1.02 because of rounding errors in the TRANSFAC data
        //    where it is only given to 2 sig fig
        double total = 0;
        BOOST_FOREACH( nucleo_dist & c, counts )
        {
            double count = c.get_total();
            total += count;
            max_count  = max(max_count,count);
            min_count  = min(min_count,count);
        }
    

        if (n_sites) {
            double average_count = total/counts.size();
            //    IF the reported number of sites is out of kilter with the raw count
            //    data then adjust the pseudo count.
            if ((n_sites > (1.2 * max_count)) || (n_sites < (0.8 * min_count))) {
                multiplier = average_count/n_sites;
            }
        } else {
            if (max_count <= min_count * 1.05) {
                n_sites = max_count;
            } else if (max_count <= min_count + 1) {
                n_sites = max_count;
            }

            //Look to see if the data is non integer implying that we have the equivalent of data from
            //lots of sites.  If so adjust the pseudo count on the assumption that it is for 100
            //sites

            for (int i = 0; (i < int(counts.size())) && (multiplier == 1); i++ ) {

                const nucleo_dist & c = counts[i];
                for (int j = 0;j < 4;j++) {

                    double frac = fmod(c.get(j),1.0);
                    if ((frac > 0.1) && (frac < 0.9)) {
                        multiplier = 1.0f/ANALOGUE_SITE_EQUIVALENT;
                        n_sites = ANALOGUE_SITE_EQUIVALENT;
                        break;
                    }
                }
            }
        }

        double pseudo_count = pssm_parameters::singleton().pseudo_counts * multiplier;
        nucleo_dist::vec dists( counts + pseudo_count );
        pssm_ptr _pssm = create_pssm( dists );

        return pssm_info(
            counts,
            pseudo_count,
            n_sites,
            _pssm,
            calculate_likelihoods_under_pssm( _pssm, dists ),
            calculate_likelihoods_under_background( _pssm )
        );
    }
};


pssm_type
get_pssm_type( const std::string & id ) {
    if( regex_match( id, PssmFromName::transfac_re ) ) return PSSM_TRANSFAC;
    if( regex_match( id, PssmFromName::custom_re ) ) return PSSM_CUSTOM;
    return PSSM_UNKNOWN_TYPE;
}

const regex PssmFromName::transfac_re( "[MR][0-9][0-9][0-9][0-9][0-9]" );
const regex PssmFromName::custom_re( "([A-Z]+)-([0-9]+)" );

boost::filesystem::path
get_pssm_cache_serialised_file()
{
    using namespace boost::filesystem;
    return
        path( BIO_NS::BioEnvironment::singleton().get_serialised_dir() )
        / path( "pssm_cache.txt" );
}

/**
 * Our cache of PSSM information
 */
struct pssm_cache
    : BIO_NS::Cache< PssmFromName >
    , BIO_NS::Singleton< pssm_cache >
{
    void init_singleton()
    {
        BIO_NS::try_to_deserialise< false >( this->elements, get_pssm_cache_serialised_file() );

        //    We don't persist the dists because they are derivable from the counts
        BOOST_FOREACH(map_t::value_type & t,elements)
        {
            t.second._dists = t.second._counts + t.second._pseudo_count;
        }
    }

    void serialise() const
    {
        BIO_NS::serialise< false >( this->elements, get_pssm_cache_serialised_file() );
    }
};

} //namespace detail




likelihoods_ptr
calculate_likelihoods_under_pssm( pssm_ptr _pssm, const nucleo_dist::vec & dists ) {
    likelihoods_ptr under_pssm( new likelihoods( pssm_parameters::singleton().likelihoods_size ) );
    calculate_pssm_likelihoods( *_pssm, dists, *under_pssm, pssm_parameters::singleton().calculate_likelihoods_map_size, false );
    //std::cout << "Under pssm:\n";
    //BOOST_FOREACH( double p, *under_pssm )
    //{
    //    std::cout << p << "\n";
    //}
    return under_pssm;
}


likelihoods_ptr
calculate_likelihoods_under_background( pssm_ptr _pssm ) {
    likelihoods_ptr under_background( new likelihoods( pssm_parameters::singleton().likelihoods_size ) );
    const nucleo_dist::vec background_dist =
        boost::assign::list_of( uniform_nucleo_dist() ).repeat( _pssm->size() - 1, uniform_nucleo_dist() );
    calculate_pssm_likelihoods( *_pssm, background_dist, *under_background, pssm_parameters::singleton().calculate_likelihoods_map_size, false );
    //std::cout << "Under background:\n";
    //BOOST_FOREACH( double p, *under_background )
    //{
    //    std::cout << p << "\n";
    //}
    return under_background;
}



bool is_transfac_pssm( const std::string & pssm_name )
{
    return detail::PSSM_TRANSFAC == detail::get_pssm_type( pssm_name );
}



bool
add_pssm_to_cache( const std::string & name, const pssm_info & pssm ) {
    return detail::pssm_cache::singleton().insert( name, pssm );
}

void
save_pssm_cache_state( )
{
    detail::pssm_cache::singleton().serialise();
}

void
clear_pssm_cache( )
{
    detail::pssm_cache::singleton().clear();
}

const pssm_info &
get_pssm(
    const std::string & name )
{
    return detail::pssm_cache::singleton()( name );
}

std::string get_pssm_name( const std::string & id )
{
    using namespace detail;
    switch( get_pssm_type( id ) )
    {
    case PSSM_TRANSFAC: return BIO_NS::parse_table_link_accession_number( id ).get_name();
    case PSSM_CUSTOM: return get_custom_pssm( id )->name;
    default: throw std::logic_error( BIOPSY_MAKE_STRING( "Do not know about PSSM with id: " << id ) );
    }
}

std::string get_pssm_url( const std::string & id )
{
    using namespace detail;
    switch( get_pssm_type( id ) )
    {
    case PSSM_TRANSFAC: return BIO_NS::parse_table_link_accession_number( id ).get_url();
    case PSSM_CUSTOM: return get_custom_pssm( id )->url;
    default: throw std::logic_error( BIOPSY_MAKE_STRING( "Do not know about PSSM with id: " << id ) );
    }
}


typedef std::vector< double > likelihoods;

unsigned get_likelihood_index( unsigned size, double score )
{
    if (size == 0)
    {
        throw std::logic_error( "size == 0" );
    }

    if (0.0 > score)
    {
        throw std::logic_error( BIOPSY_MAKE_STRING( "Score < 0: " << score ) );
    }

    if (score * 2 * size > 2 * size + 1)
    {
        throw std::logic_error( BIOPSY_MAKE_STRING( "Score > 1 and a bit: " << score ) );
    }

    unsigned idx = ( unsigned )(score * size);
    if( size == idx ) //cater for a perfect match
    {
        --idx;
    }
    BOOST_ASSERT( 0 <= idx && idx < size );

    return idx;
}

double get_likelihood( likelihoods_ptr likelihoods, double score )
{
    const unsigned index = get_likelihood_index( likelihoods->size(), score );
    return ( *likelihoods )[ index ];
}


likelihoods_ptr
accumulate_likelihoods( likelihoods_ptr ls )
{
    const unsigned size = ls->size();
    likelihoods_ptr result( new likelihoods( size ) );
    double cum = 0.0;
    for( unsigned i = 0; size != i; ++i )
    {
        const unsigned idx = size - 1 - i;
        cum += ( *ls )[ idx ];
        ( *result )[ idx ] = cum;
    }
    return result;
}


double
normalise_likelihoods( likelihoods & result )
{
    const double sum = std::accumulate( result.begin(), result.end(), 0.0 );
    if( 0 >= sum )
    {
        throw std::logic_error( BIOPSY_MAKE_STRING( "Cannot normalise likelihood vector with sum < 0: " << sum ) );
    }

    BOOST_FOREACH( double & p, result )
    {
        p /= sum;
    }

    return sum;
}


void
calculate_pssm_likelihoods(
    const pssm & pssm,
    const nucleo_dist::vec & distribution,
    likelihoods & result,
    const size_t max_map_size,
    bool verbose )
{
    using namespace boost::numeric;

    /** Maps scores to probabilities. */
    typedef std::map< double, double > prob_map;
    typedef boost::scoped_ptr< prob_map > prob_map_ptr;

    if( pssm.size() != distribution.size() )
    {
        BOOST_FOREACH( const nucleo_dist & d, distribution )
        {
            std::cout << d.get( 0 ) << "," << d.get( 1 ) << "," << d.get( 2 ) << "," << d.get( 3 ) << "\n";
        }

        throw
            std::invalid_argument(
                BIOPSY_MAKE_STRING(
                    "pssm and distribution are different sizes: "
                    << pssm.size() << " != " << distribution.size() ) );
    }

    if( result.empty() )
    {
        throw std::invalid_argument( "result vector is empty" );
    }

    prob_map_ptr probs1( new prob_map );
    prob_map_ptr probs2( new prob_map );
    prob_map * last_probs = probs1.get(); //these alternate
    prob_map * new_probs = probs2.get();
    (*last_probs)[0.0] = std::log(1.0); //initialise with score of 0 has likelihood 1.0 at the start of the algorithm

    //for each row
    for( unsigned base = 0; pssm.size() != base; ++base )
    {
        new_probs->clear(); //start afresh

        //go to the next nucleotide if no observations
        if ( 0 == distribution[ base ].get_total() )
        {
            continue;
        }

        //if we have too many scores we need to quantise them
        while (last_probs->size() > max_map_size)
        {
            if( verbose )
            {
                std::cout << "Reducing map size from " <<  last_probs->size();
            }

            const double min_score = last_probs->begin()->first;
            const double max_score = last_probs->rbegin()->first;
            const double quantum = (max_score - min_score) / max_map_size;

            //for each pair of consecutive scores
            for( prob_map::iterator p1 = last_probs->begin(); last_probs->end() != p1; )
            {
                //get the next score
                prob_map::iterator p2 = p1;
                ++p2;
                if( last_probs->end() == p2 )
                {
                    break;
                }

                const double score1 = p1->first;
                const double score2 = p2->first;
                BOOST_ASSERT(score2 > score1); //problem if in wrong order...

                //are the two scores close?
                if( score2 - score1 < quantum )
                {
                    //we remove them both and insert a replacement for the average score that is as likely
                    //as both together
                    prob_map::iterator to_erase1 = p1;
                    prob_map::iterator to_erase2 = p2;

                    //insert the new value
                    const double prob1 = std::exp( p1->second );
                    const double prob2 = std::exp( p2->second );
                    const double prob_sum = prob1 + prob2;
#ifndef NDEBUG
                    const double relative_prob1 = prob1 / prob_sum;
#endif
                    const double relative_prob2 = prob2 / prob_sum;
                    const double avg_score = score1 + ( score2 - score1 ) * relative_prob2; //had to be careful here with floating point underflow
                    BOOST_ASSERT( score1 <= avg_score );
                    BOOST_ASSERT( avg_score <= score2 );

                    //move on before we insert the new value and erase the old
                    p1 = p2; ++p1;

                    //erase the old values
                    last_probs->erase(to_erase1);
                    last_probs->erase(to_erase2);

                    //insert the new value
                    const double log_prob_sum = std::log( prob_sum );
                    std::pair< prob_map::iterator, bool > insert_result =
                        last_probs->insert( prob_map::value_type( avg_score, log_prob_sum ) );

                    //did it actually insert?
                    if( ! insert_result.second )
                    {
                        //no - so update the existing element
                        insert_result.first->second =
                            std::log( std::exp( insert_result.first->second ) + std::exp( log_prob_sum ) );
                    }
                }
                else
                {
                    //move on
                    p1 = p2;
                }
            }

            if (verbose)
            {
                std::cout << " to " <<  last_probs->size() << std::endl;
            }

        } //while map is too large

        //for each nucleotide
        for( unsigned i = 0; 4 != i; ++i )
        {
            //likelihood in the pssm for this nucleotide in this position? Use a psuedo count of 1
            const double prob_this_nucleo = distribution[ base ].get_freq( i );

            //if none, cannot contribute to the total probabilities
            if( 0.0 == prob_this_nucleo )
            {
                continue;
            }
            //BOOST_ASSERT(0 != entry->get_count(*n)); - this was only true before pseudo counts used

            const double log_prob_this_nucleo = std::log( prob_this_nucleo );
            const double score_this_nucleo = pssm[ base ].get( i );

            //for each score we have already achieved
            for( prob_map::const_iterator p = last_probs->begin(); last_probs->end() != p; ++p )
            {
                const double new_score = p->first + score_this_nucleo;
                const double new_log_likelihood = p->second + log_prob_this_nucleo;

                //update the new probability map

                //first try to insert
                std::pair< prob_map::iterator, bool > insert_result =
                    new_probs->insert( prob_map::value_type( new_score, new_log_likelihood ) );

                //did it actually insert?
                if( ! insert_result.second )
                {
                    //no - so update the existing element
                    insert_result.first->second =
                        std::log( std::exp( insert_result.first->second ) + std::exp( new_log_likelihood ) );
                }
            }

        } //for each nucleotide

        if( verbose )
        {
            //for debugging
            double total_prob = 0.0;
            for( prob_map::const_iterator p = new_probs->begin(); new_probs->end() != p; ++p )
            {
                total_prob += std::exp( p->second );
            }
            std::cout
                << "Base " << base << " has prob map size " << new_probs->size()
                << " and total probability: " << total_prob
                << std::endl;

        }

        std::swap( new_probs, last_probs );

    } //for each row

    //initialise the result vector
    std::fill( result.begin(), result.end(), 0.0 );

    //for each score in the last_probs map, place in the likelihoods vector
    for( prob_map::const_iterator p = last_probs->begin(); last_probs->end() != p; ++p )
    {
        //normalise the score
        BOOST_ASSERT( in( p->first, interval< double >( -0.0001, 1.0001 ) ) );
        const double score = std::min( 1.0, std::max( 0.0, p->first ) );

        //what is the probability
        const double log_likelihood = p->second;
        const double prob = std::exp( log_likelihood );

        result[ get_likelihood_index( result.size(), score ) ] += prob;
    }

    normalise_likelihoods( result );
}

double
get_p_binding_using_p_value(
    double score,
    likelihoods_ptr likelihoods )
{
    const unsigned idx = get_likelihood_index( likelihoods->size(), score );
    const double likelihood_under_background = ( *likelihoods )[ idx ];
    const double likelihood_under_binding = 1.0 / double( likelihoods->size() );
    return
        get_p_binding(
            get_odds_ratio(
                likelihood_under_binding,
                likelihood_under_background ) );
}

double
get_p_binding_from_score(
    const pssm_info & p,
    double score )
{
    const pssm_parameters & params = pssm_parameters::singleton();
    const bool use_cumulative_dists = params.use_cumulative_dists;
    if( params.use_p_value )
    {
        return
            get_p_binding_using_p_value(
                score,
                p.get_dist( false, use_cumulative_dists ) );
    }
    else
    {
        return
            get_p_binding(
                get_odds_ratio(
                    score,
                    p.get_dist( true, use_cumulative_dists ),
                    p.get_dist( false, use_cumulative_dists ) ) );
    }
}


const pssm_info::matrix_t &
pssm_info::get_log_likelihoods() const {

    // check whether we have calculated them
    if( ! _log_likelihoods ) {
        _log_likelihoods.reset( new matrix_t( boost::extents[ boost::size( _dists ) ][ 4 ] ) );
        for( size_t i = 0; _dists.size() != i; ++i ) {
            for( size_t b = 0; 4 != b; ++b ) {
                ( *_log_likelihoods )[ i ][ b ] = std::log( _dists[ i ].get_freq( b ) );
            }
        }
    }

    return *_log_likelihoods;
}



double
get_p_binding_on_sequence(
    const pssm_info & p,
    sequence::const_iterator s_begin )
{
    return get_p_binding_from_score(
        p,
        score( *( p._pssm ), s_begin ) );
}


/**
Scores a pssm on the reverse complement of the sequence and adjusts for distributions.
*/
double
get_p_binding_on_reverse_complement(
    const pssm_info & p,
    sequence::const_iterator s_begin )
{
    return
        get_p_binding_from_score(
            p,
            score_complement( *( p._pssm ), s_begin ) );
}


pssm_info::pssm_info(
    const nucleo_dist::vec & counts,
    double pseudo_count,
    int number_of_sites,
    pssm_ptr p,
    likelihoods_ptr binding_dist,
    likelihoods_ptr background_dist )
    : _counts( counts )
    , _dists( counts + pseudo_count )
    , _pseudo_count( pseudo_count )
    , _number_of_sites( number_of_sites )
    , _pssm( p )
    , _binding_dist( binding_dist )
    , _background_dist( background_dist )
{
}

likelihoods_ptr
pssm_info::get_dist( bool binding, bool cumulative ) const
{
    if( cumulative )
    {
        if( binding )
        {
            if( ! _cumulative_binding_dist )
            {
                _cumulative_binding_dist = accumulate_likelihoods( get_dist( binding, false ) );
            }
            return _cumulative_binding_dist;
        }
        else
        {
            if( ! _cumulative_background_dist )
            {
                _cumulative_background_dist = accumulate_likelihoods( get_dist( binding, false ) );
            }
            return _cumulative_background_dist;
        }
    }
    else
    {
        if( binding )
        {
            if( ! _binding_dist )
            {
                throw std::logic_error( "Do not have binding distribution." );
            }
            return _binding_dist;
        }
        else
        {
            if( ! _background_dist )
            {
                throw std::logic_error( "Do not have background distribution." );
            }
            return _background_dist;
        }
    }
}


/** Convert a biopsy pssm to a bio pssm. */
BIO_NS::Pssm make_transfac_pssm( const pssm & _pssm )
{
    BIO_NS::Pssm result;
    BOOST_FOREACH( const nucleo_dist & dist, _pssm )
        result.push_back( BIO_NS::PssmEntry( (BIO_NS::float_t)(10*dist.get(0)),(BIO_NS::float_t)(10*dist.get(1)), 
         (BIO_NS::float_t)(10*dist.get(2)), (BIO_NS::float_t)(10*dist.get(3))) );
    return result;
}




} //namespace biopsy
