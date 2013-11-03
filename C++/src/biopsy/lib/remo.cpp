/**
@file

Copyright John Reid 2006

*/

#include "biopsy/defs.h"
#include "biopsy/remo.h"

#include <bio/remo.h>
#include <bio/serialisable.h>


namespace boost {
namespace serialization {

template< typename Archive >
void serialize(
    Archive & ar,
    biopsy::remo::location & v,
    const unsigned int version )
{
    ar & v._start;
    ar & v._end;
}

template< typename Archive >
void serialize(
    Archive & ar,
    biopsy::remo::ensembl_id & v,
    const unsigned int version )
{
    ar & v._prefix;
    ar & v._num;
}

template< typename Archive >
void serialize(
    Archive & ar,
    biopsy::remo::exon & v,
    const unsigned int version )
{
    ar & v._location;
    ar & v._id;
}

template< typename Archive >
void serialize(
    Archive & ar,
    biopsy::remo::ensembl_database_id & v,
    const unsigned int version )
{
    ar & v._species;
    ar & v._software_version;
    ar & v._ncbi_build;
    ar & v._build_version;
}

template< typename Archive >
void serialize(
    Archive & ar,
    biopsy::remo::alignment_sequence_id & v,
    const unsigned int version )
{
    ar & v._gene_id;
    ar & v._transcript_id;
    ar & v._db_id;
    ar & v._version;
}

template< typename Archive >
void serialize(
    Archive & ar,
    biopsy::remo::alignment_sequence_info & v,
    const unsigned int version )
{
    ar & v._length;
    ar & v._has_position;
    ar & v._position;
    ar & v._exons;
    ar & v._region;
}

template< typename Archive >
void serialize(
    Archive & ar,
    biopsy::remo::remo_sequence & v,
    const unsigned int version )
{
    ar & v._masked_sequence;
    ar & v._unmasked_sequence;
    ar & v._location;
    ar & v._target_location;
    ar & v._conservation;
    ar & v._repeat_ratio;
    ar & v._belief;
}

template< typename Archive >
void serialize(
    Archive & ar,
    biopsy::remo::module & v,
    const unsigned int version )
{
    ar & v._sequences;
}

template< typename Archive >
void serialize(
    Archive & ar,
    biopsy::remo::aligned_sequence_set & v,
    const unsigned int version )
{
    ar & v._sequences;
    ar & v._centre_sequence;
}

template< typename Archive >
void serialize(
    Archive & ar,
    biopsy::remo::ensembl_id_alignment_map & v,
    const unsigned int version )
{
    ar & v._map;
}

template< typename Archive >
void serialize(
    Archive & ar,
    biopsy::remo::remome & v,
    const unsigned int version )
{
    ar & v._remos;
}

} // namespace serialization
} // namespace boost




namespace biopsy {
namespace remo {

location::location( )
: _start( -1 )
, _end( -1 )
{
}

location::location( int start, int end )
: _start( start )
, _end( end )
{
}

bool location::operator==( const location & rhs ) const
{
    return
        ( _start == rhs._start )
        &&
        ( _end == rhs._end )
        ;
}

bool location::operator<( const location & rhs ) const
{
    if( _start < rhs._start ) return true;
    if( rhs._start < _start ) return false;
    if( _end < rhs._end ) return true;
    return false;
}

std::string
location::str() const
{
    return BIOPSY_MAKE_STRING( *this );
}



std::ostream &
operator<<( std::ostream & os, const location & l )
{
    os << '[' << l._start << ',' << l._end << ']';
    return os;
}


std::ostream &
operator<<( std::ostream & os, region r )
{
    switch( r )
    {
    case region_upstream: return os << "upstream";
    case region_downstream: return os << "downstream";
    case region_gene: return os << "gene region";
    case region_undefined: return os << "<undefined region>";
    default: throw std::invalid_argument( "Unknown region enum" );
    }
}

ensembl_id::ensembl_id( const std::string & prefix, unsigned num ) : _prefix( prefix ), _num( num ) { }

std::string
ensembl_id::str() const
{
    return BIOPSY_MAKE_STRING( *this );
}


bool ensembl_id::looks_like( const std::string & text )
{
    static const boost::regex re( "[A-Z]+[0-9]+" );
    return boost::regex_match( text.c_str(), re );
}

ensembl_id ensembl_id::parse( const std::string & text )
{
    static const boost::regex re( "([A-Z]+)([0-9]+)" );
    boost::cmatch what;
    const bool result = boost::regex_match( text.c_str(), what, re );
    if( result )
    {
        return
            ensembl_id(
                what[ 1 ],
                boost::lexical_cast< unsigned >( what[ 2 ] ) );
    }
    else
    {
        throw std::logic_error( BIOPSY_MAKE_STRING( "Could not parse ensembl_id: \"" << text << "\"" ) );
    }
}

bool ensembl_id::operator==( const ensembl_id & rhs ) const
{
    return _prefix == rhs._prefix && _num == rhs._num;
}

bool ensembl_id::operator<( const ensembl_id & rhs ) const
{
    if( _prefix < rhs._prefix ) return true;
    if( rhs._prefix < _prefix ) return false;
    if( _num < rhs._num ) return true;
    return false;
}

std::ostream &
operator<<( std::ostream & os, const ensembl_id & id )
{
    boost::io::ios_fill_saver  ios_saver( os );
    os.fill( '0' );

    return os
        << id._prefix
        << std::setw( 11 ) << id._num
        ;
}

exon::exon( )
{
}

exon::exon( location l, const ensembl_id & id ) : _location( l ), _id( id ) { }

bool exon::operator==( const exon & rhs ) const
{
    return _location == rhs._location && _id == rhs._id;
}

bool exon::operator<( const exon & rhs ) const
{
    if( _location < rhs._location ) return true;
    if( rhs._location < _location ) return false;
    if( _id < rhs._id ) return true;
    return false;
}

ensembl_database_id::ensembl_database_id( )
: _software_version( -1 )
, _ncbi_build( -1 )
{
}

ensembl_database_id::ensembl_database_id(
    const species & s,
    int software_version,
    int ncbi_build,
    const std::string & build_version )
    : _species( s )
    , _software_version( software_version )
    , _ncbi_build( ncbi_build )
    , _build_version( build_version )
{
}

ensembl_database_id::ensembl_database_id( const std::string & id )
    : _species( "" )
    , _software_version( 0 )
    , _ncbi_build( 0 )
    , _build_version( id )
{
}

bool
ensembl_database_id::operator<( const ensembl_database_id & rhs ) const
{
    if( _species < rhs._species ) return true;
    if( rhs._species < _species ) return false;
    if( _software_version < rhs._software_version ) return true;
    if( rhs._software_version < _software_version ) return false;
    if( _ncbi_build < rhs._ncbi_build ) return true;
    if( rhs._ncbi_build < _ncbi_build ) return false;
    if( _build_version < rhs._build_version ) return true;
    return false;
}


bool
ensembl_database_id::operator==( const ensembl_database_id & rhs ) const
{
    return
        ( _species == rhs._species )
        && ( _software_version == rhs._software_version )
        && ( _ncbi_build == rhs._ncbi_build )
        && ( _build_version == rhs._build_version )
        ;
}

std::string
ensembl_database_id::str() const
{
    return BIOPSY_MAKE_STRING( *this );
}



namespace detail
{
    //want to match this sort of thing
    //        mus_musculus_core_28_33d
    static boost::regex ensembl_database_id_re( "([a-z_A-Z]+)_core_([0-9]+)_([0-9]+)([a-zA-Z]*)" );
}

bool
ensembl_database_id::looks_like( const std::string & tag )
{
    return boost::regex_match( tag.c_str(), detail::ensembl_database_id_re );
}

ensembl_database_id
ensembl_database_id::parse( const std::string & tag )
{
    boost::cmatch what;
    const bool result = boost::regex_match( tag.c_str(), what, detail::ensembl_database_id_re );
    if( result )
    {
        return
            ensembl_database_id(
                what[ 1 ],
                boost::lexical_cast< int >( what[ 2 ] ),
                boost::lexical_cast< int >( what[ 3 ] ),
                what[ 4 ] );
    }
    else
    {
        throw std::invalid_argument( BIOPSY_MAKE_STRING( "Could not parse: " << tag ) );
    }
}



std::ostream &
operator<<( std::ostream & os, const ensembl_database_id & id )
{
    os
        << id._species
        << "_core_" << id._software_version
        << '_' << id._ncbi_build
        << id._build_version
        ;

    return os;
}

alignment_sequence_id::alignment_sequence_id()
: _version( -1 )
{
}

alignment_sequence_id::alignment_sequence_id( const std::string & id )
: _db_id( id )
, _version( 0 )
{
}


alignment_sequence_id::alignment_sequence_id(
    const ensembl_id & gene_id,
    const ensembl_id & transcript_id,
    const ensembl_database_id & db_id,
    int version )
    : _gene_id( gene_id )
    , _transcript_id( transcript_id )
    , _db_id( db_id )
    , _version( version )
{
}


bool alignment_sequence_id::operator==( const alignment_sequence_id & rhs ) const
{
    return
        ( _gene_id == rhs._gene_id )
        && ( _transcript_id == rhs._transcript_id )
        && ( _db_id == rhs._db_id )
        && ( _version == rhs._version )
        ;
}

bool alignment_sequence_id::operator<( const alignment_sequence_id & rhs ) const
{
    if( _gene_id < rhs._gene_id ) return true;
    if( rhs._gene_id < _gene_id ) return false;
    if( _transcript_id < rhs._transcript_id ) return true;
    if( rhs._transcript_id < _transcript_id ) return false;
    if( _db_id < rhs._db_id ) return true;
    if( rhs._db_id < _db_id ) return false;
    if( _version < rhs._version ) return true;
    return false;
}

std::string
alignment_sequence_id::str() const
{
    return BIOPSY_MAKE_STRING( *this );
}

namespace detail
{
    static boost::regex alignment_sequence_id_re( "([A-Z0-9]+) ([A-Z0-9]+) \\(([a-zA-Z_0-9]+)\\)(?: ([0-9]+))?" );
}

bool
alignment_sequence_id::looks_like( const std::string & tag )
{
    return boost::regex_match( tag.c_str(), detail::alignment_sequence_id_re );
}

alignment_sequence_id::ptr
alignment_sequence_id::parse( const std::string & tag )
{
    boost::cmatch what;
    const bool result = boost::regex_match( tag.c_str(), what, detail::alignment_sequence_id_re );
    if( result )
    {
        return
            alignment_sequence_id::ptr(
                new alignment_sequence_id(
                    ensembl_id::parse( what[ 1 ] ),
                    ensembl_id::parse( what[ 2 ] ),
                    ensembl_database_id::parse( what[ 3 ] ),
                    what[4].matched
                        ? boost::lexical_cast< int >( what[ 4 ] )
                        : 0 ) );
    }
    else
    {
        throw std::invalid_argument( BIOPSY_MAKE_STRING( "Could not parse: " << tag ) );
    }
}



std::ostream &
operator<<( std::ostream & os, const alignment_sequence_id & id )
{
    os
        << id._gene_id
        << ' ' << id._transcript_id
        << " (" << id._db_id
        << ")"
        ;

    if ( 0 != id._version )
    {
        os << " " << id._version;
    }

    return os;
}

alignment_sequence_info::alignment_sequence_info( )
: _length( 0 )
, _has_position( false )
, _position( 0 )
, _region( region_undefined )
{
}


alignment_sequence_info::alignment_sequence_info(
    int length,
    bool has_position,
    int position,
    exon::list_ptr exons,
    region r )
    : _length( length )
    , _has_position( has_position )
    , _position( position )
    , _exons( exons )
    , _region( r )
{
}

remo_sequence::remo_sequence( )
: _conservation( 0 )
, _repeat_ratio( 0 )
, _belief( 0 )
{
}


remo_sequence::remo_sequence(
    const sequence & masked_sequence,
    const sequence & unmasked_sequence,
    location l,
    location target_l,
    unsigned conservation,
    unsigned repeat_ratio,
    double belief )
    : _masked_sequence( masked_sequence )
    , _unmasked_sequence( unmasked_sequence )
    , _location( l )
    , _target_location( target_l )
    , _conservation( conservation )
    , _repeat_ratio( repeat_ratio )
    , _belief( belief )
{
}

const sequence &
remo_sequence::get_sequence( bool masked ) const
{
    return masked ? _masked_sequence : _unmasked_sequence;
}

bool
remo_sequence::operator==( const remo_sequence & rhs) const
{
    if( _location != rhs._location ) return false;
    if( _target_location != rhs._target_location ) return false;
    if( _conservation != rhs._conservation ) return false;
    if( _repeat_ratio != rhs._repeat_ratio ) return false;
    if( _masked_sequence != rhs._masked_sequence ) return false;
    if( _unmasked_sequence != rhs._unmasked_sequence ) return false;
    if( _belief != rhs._belief ) return false;
    return true;
}

bool
remo_sequence::operator<( const remo_sequence & rhs) const
{
    if( _location < rhs._location ) return true;
    if( rhs._location < _location ) return false;
    if( _target_location < rhs._target_location ) return true;
    if( rhs._target_location < _target_location ) return false;
    if( _conservation < rhs._conservation ) return true;
    if( rhs._conservation < _conservation ) return false;
    if( _repeat_ratio < rhs._repeat_ratio ) return true;
    if( rhs._repeat_ratio < _repeat_ratio ) return false;
    if( _masked_sequence < rhs._masked_sequence ) return true;
    if( rhs._masked_sequence < _masked_sequence ) return false;
    if( _unmasked_sequence < rhs._unmasked_sequence ) return true;
    if( rhs._unmasked_sequence < _unmasked_sequence ) return false;
    if( _belief < rhs._belief ) return true;
    return false;
}

bool
module::operator==( const module & rhs ) const
{
    if( _sequences != _sequences ) return false;
    return true;
}

bool
module::operator<( const module & rhs ) const
{
    if( _sequences < _sequences ) return true;
    return false;
}

alignment_sequence_id::list_ptr
module::get_sequence_ids() const
{
    alignment_sequence_id::list_ptr result( new alignment_sequence_id::list );
    BOOST_FOREACH( sequence_remo_map::value_type v, _sequences )
    {
        result->push_back( v.first );
    }
    return result;
}

remo_sequence::list_ptr
module::get_sequences( alignment_sequence_id::ptr id ) const
{
    sequence_remo_map::const_iterator i = _sequences.find( id );
    if( _sequences.end() == i )
    {
        throw std::logic_error( "Could not find sequence in remo" );
    }
    return i->second;
}

sequence
module::get_sequence_for( alignment_sequence_id::ptr id, bool masked ) const
{
    sequence_remo_map::const_iterator i = _sequences.find( id );
    if( _sequences.end() == i )
    {
        throw std::logic_error( "Could not find sequence in remo" );
    }
    sequence result;
    //int end_of_previous_sequence = std::numeric_limits< int >::min();
    BOOST_FOREACH( remo_sequence::ptr s, *( i->second ) )
    {
        //BOOST_ASSERT( s->_location._start > end_of_previous_sequence );
        //end_of_previous_sequence = s->_location._end;

        result += s->get_sequence( masked );
    }
    return result;
}


alignment_sequence_id::list_ptr
aligned_sequence_set::get_sequence_ids() const
{
    alignment_sequence_id::list_ptr result( new alignment_sequence_id::list );
    BOOST_FOREACH( const alignment_sequence_info_map::value_type & v, _sequences )
    {
        result->push_back( v.first );
    }
    return result;
}

alignment_sequence_info::ptr
aligned_sequence_set::get_sequence_info( alignment_sequence_id::ptr seq_id ) const
{
    alignment_sequence_info_map::const_iterator i = _sequences.find( seq_id );
    if( _sequences.end() == i )
    {
        throw std::logic_error( "Cannot find sequence id" );
    }
    return i->second;
}

bool
aligned_sequence_set::operator==( const aligned_sequence_set & rhs ) const
{
    return *( _centre_sequence ) == *( rhs._centre_sequence ) && _sequences == rhs._sequences;
}

bool
aligned_sequence_set::operator<( const aligned_sequence_set & rhs ) const
{
    if( *( _centre_sequence ) < *( rhs._centre_sequence ) ) return true;
    if( *( rhs._centre_sequence ) < *( _centre_sequence ) ) return false;
    if( _sequences < rhs._sequences ) return true;
    return false;
}

ensembl_id::list_ptr
ensembl_id_alignment_map::get_genes( const std::string & prefix ) const
{
    ensembl_id::list_ptr result( new ensembl_id::list );
    BOOST_FOREACH( const ensembl_aligned_sequences_map::value_type & v, _map )
    {
        if( "" == prefix || prefix == v.first._prefix )
        {
            result->push_back( v.first );
        }
    }
    return result;
}

aligned_sequence_set::list_ptr
ensembl_id_alignment_map::get_alignments_for( const ensembl_id & id ) const
{
    ensembl_aligned_sequences_map::const_iterator i = _map.find( id );
    return
        _map.end() == i
            ? aligned_sequence_set::list_ptr( new aligned_sequence_set::list )
            : i->second;
}

namespace detail {


alignment_sequence_id::ptr
alignment_sequence_id_from_extraction_sequence_tag( const std::string & tag )
{
    USING_BIO_NS;

    alignment_sequence_id::ptr result;

    if( alignment_sequence_id::looks_like( tag ) )
    {
        return alignment_sequence_id::parse( tag );
    }
    else
    {
        result.reset(
            new alignment_sequence_id(
                tag ) );
    }

    return result;
}

location
get_location( const BIO_NS::ReMoRange & range )
{
    return location( range.start, range.end );
}

ensembl_id
get_ensembl_id( const BIO_NS::EnsemblId & id )
{
#if 0 //have fixed bug..
    //there is a bug in the BIO_NS parsing so check if the value has a number in it, if so reparse it
    static boost::regex re( "[0-9]" );
    if( boost::regex_search( id.value.c_str(), re ) )
    {
        return ensembl_id::parse( BIOPSY_MAKE_STRING( id.value << id.number ) );
    }
#endif
    return ensembl_id( id.value, id.number );
}

exon::list_ptr
get_exons( const BIO_NS::ReMoExon::list_t & exons )
{
    exon::list_ptr result( new exon::list );
    BOOST_FOREACH( const BIO_NS::ReMoExon & e, exons )
    {
        result->push_back( exon( get_location( e.range ), get_ensembl_id( e.id ) ) );
    }
    return result;
}

region
get_region( BIO_NS::ReMoLocation location )
{
    switch( location )
    {
    case BIO_NS::REMO_LOC_UPSTREAM: return region_upstream;
    case BIO_NS::REMO_LOC_DOWNSTREAM: return region_downstream;
    case BIO_NS::REMO_LOC_GENEREGION: return region_gene;
    case BIO_NS::REMO_LOC_UNDEFINED: return region_undefined;
    default:
        throw std::invalid_argument( "Unknown ReMoLocation" );
    }
}

template< typename T, typename Cmp >
struct ensure_unique
{
    std::set< T, Cmp > _universe;
    T operator()( T t )
    {
        return *( _universe.insert( t ).first );
    }
};

/** make sure we don't have duplicate sequence ids hanging around */
static ensure_unique< alignment_sequence_id::ptr, smart_ptr_less_than_value > _ensure_unique;

remome::ptr
create_remome_from_bio_remo_extraction( BIO_NS::ReMoExtraction::ptr_t extraction )
{
    USING_BIO_NS;

    remome::remo_map_ptr remo_map( new remome::remo_map );
    BOOST_FOREACH( ReMoSequenceGroup::ptr_t group, extraction->sequence_groups )
    {
        //get the sequences that make up this aligned group
        aligned_sequence_set::ptr _aligned_sequences( new aligned_sequence_set );
        BOOST_FOREACH( ReMoSequence::ptr_t seq, group->sequences )
        {
#if 0
            std::cout
                << seq->id << "\n"
                << seq->length << "\n"
                << seq->location << "\n"
                << seq->has_position << "\n"
                << seq->position << "\n"
                << seq->species << "\n"
                << "\n";
#endif
            //create a sequence id from the details
            alignment_sequence_id::ptr seq_id = _ensure_unique( alignment_sequence_id_from_extraction_sequence_tag( seq->id ) );

            //create the info from the details
            alignment_sequence_info::ptr info(
                new alignment_sequence_info(
                    seq->length,
                    seq->has_position,
                    seq->position,
                    get_exons( seq->exons ),
                    get_region( seq->location) ) );

            //put it in the map
#ifndef NDEBUG
            const bool was_inserted =
#endif
            _aligned_sequences->_sequences.insert(
                aligned_sequence_set::alignment_sequence_info_map::value_type(
                    seq_id,
                    info ) ).second;
            BOOST_ASSERT( was_inserted );
        }

        //get the remos for each one
        module::list_ptr remos( new module::list );
        BOOST_FOREACH( const ReMoBundle::map_t::value_type & v, group->remo_bundles )
        {
            ReMoBundle::ptr_t bundle = v.second;

            //get the centre sequence from the details
            alignment_sequence_id::ptr centre_seq = _ensure_unique( alignment_sequence_id_from_extraction_sequence_tag( bundle->centre_sequence ) );
            if( 0 == _aligned_sequences->_centre_sequence )
            {
                //assign as the centre sequence for our sequences
                _aligned_sequences->_centre_sequence = centre_seq;
            }
            else
            {
                //check it is the same
                BOOST_ASSERT( _aligned_sequences->_centre_sequence == centre_seq );
            }


            //go through each sequence and create remo
            module::ptr remo( new module );
            BOOST_FOREACH( const ReMo::map_t::value_type & r, bundle->remos )
            {
                const std::string & seq_id = r.first;
                const ReMo::list_t & remos = r.second;

                alignment_sequence_id::ptr id = _ensure_unique( alignment_sequence_id_from_extraction_sequence_tag( seq_id ) );

                remo_sequence::list_ptr seq_remos( new remo_sequence::list );
                BOOST_FOREACH( ReMo::ptr_t s, remos )
                {
                    seq_remos->push_back(
                        remo_sequence::ptr(
                            new remo_sequence(
                                s->masked_sequence,
                                s->unmasked_sequence,
                                get_location( s->range ),
                                get_location( s->target_range ),
                                s->conservation,
                                s->repeat_ratio,
                                s->belief ) ) );
                }

#ifndef NDEBUG
                const bool was_inserted =
#endif
                remo->_sequences.insert( module::sequence_remo_map::value_type( id, seq_remos ) ).second;
                BOOST_ASSERT( was_inserted );
            }
            remos->push_back( remo );
        }
        remo_map->insert( remome::remo_map::value_type( _aligned_sequences, remos ) );
    }

    remome::ptr result( new remome );
    result->_remos = remo_map;

    return result;
}

} //namespace detail

aligned_sequence_set::list_ptr
remome::get_aligned_sequences() const
{
    aligned_sequence_set::list_ptr result( new aligned_sequence_set::list );
    BOOST_FOREACH( remo_map::value_type v, *_remos )
    {
        result->push_back( v.first );
    }
    return result;
}

module::list_ptr
remome::get_remos_for( aligned_sequence_set::ptr al ) const
{
    remo_map::const_iterator i = _remos->find( al );
    if( _remos->end() == i )
    {
        throw std::invalid_argument( "Cannot find entry for aligned sequences in remome" );
    }
    return i->second;
}

void
remome::remove_remos_for( aligned_sequence_set::ptr aligned_seqs )
{
    remo_map::iterator i = _remos->find( aligned_seqs );
    if( _remos->end() == i )
    {
        throw std::invalid_argument( "Cannot find entry for aligned sequences in remome" );
    }
    _remos->erase( i );
}

remome::ptr
remome::deserialise( const std::string & path )
{
    return
        BIO_NS::deserialise< true, remome >(
            boost::filesystem::path( path ) );
}


void
remome::serialise( const std::string & path ) const
{
    return
        BIO_NS::serialise< true >(
            *this,
            boost::filesystem::path( path ) );
}


ensembl_id_alignment_map::ptr
make_gene_alignment_map(
    aligned_sequence_set::list_ptr aligned_sequence_set )
{
    ensembl_id_alignment_map::ptr result( new ensembl_id_alignment_map );

    //for each set of aligned sequences
    BOOST_FOREACH( aligned_sequence_set::ptr aligned_seqs, *aligned_sequence_set )
    {
        //for each aligned sequence in the set
        BOOST_FOREACH( aligned_sequence_set::alignment_sequence_info_map::value_type v, aligned_seqs->_sequences )
        {
            const ensembl_id & gene = v.first->_gene_id;
            ensembl_id_alignment_map::ensembl_aligned_sequences_map::const_iterator g = result->_map.find( gene );
            if( result->_map.end() == g )
            {
                g =
                    result->_map.insert(
                        ensembl_id_alignment_map::ensembl_aligned_sequences_map::value_type(
                            gene,
                            aligned_sequence_set::list_ptr( new aligned_sequence_set::list ) ) ).first;
            }
            g->second->push_back( aligned_seqs );
        }
    }

    return result;
}

remome::ptr load_remome_from_file( const std::string & filename )
{
    USING_BIO_NS;

    const boost::filesystem::path file( filename );
    ReMoExtraction::ptr_t extraction = ReMoExtraction::deserialise( file );

    return detail::create_remome_from_bio_remo_extraction( extraction );
}


remome::ptr parse_remome_from_file( const std::string & filename )
{
    USING_BIO_NS;

    const boost::filesystem::path file( filename );
    ReMoExtraction::ptr_t extraction = parse_remo_extraction( file );

    return detail::create_remome_from_bio_remo_extraction( extraction );
}

std::string
get_remo_id( aligned_sequence_set::ptr aligned_seqs, module::ptr _remo )
{
    //check we have only one sequence in the centre sequence
    BOOST_ASSERT( 1 == _remo->get_sequences( aligned_seqs->_centre_sequence )->size() );

    return
        BIOPSY_MAKE_STRING(
            *( aligned_seqs->_centre_sequence )
            << ":" << aligned_seqs->get_sequence_info( aligned_seqs->_centre_sequence )->_region
            << ":" << ( *( _remo->get_sequences( aligned_seqs->_centre_sequence )->begin() ) )->_location )
        ;
}



remo_locator
get_remo_from_id( remome::ptr _remome, const std::string & id )
{
    static const boost::regex id_re( "([^:]+):([^:]+):([^:]+)" );
    boost::cmatch id_what;
    if( ! boost::regex_match( id.c_str(), id_what, id_re ) )
    {
        throw std::invalid_argument( BIOPSY_MAKE_STRING( "Cannot parse id: " << id ) );
    }

    //parse the location
    const std::string location_str = id_what[3];
    static const boost::regex location_re( "\\[([0-9]+),([0-9]+)\\]" );
    boost::cmatch location_what;
    if( ! boost::regex_match( location_str.c_str(), location_what, location_re ) )
    {
        throw std::invalid_argument( BIOPSY_MAKE_STRING( "Cannot parse location from id: " << id ) );
    }
    const location loc(
        boost::lexical_cast< int >( location_what[1] ),
        boost::lexical_cast< int >( location_what[2] ) );

    //get the remos for this id
    const std::string seq_tag = id_what[1];
    alignment_sequence_id::ptr seq_id = detail::_ensure_unique( detail::alignment_sequence_id_from_extraction_sequence_tag( seq_tag ) );

    //find the aligned sequences with this centre sequence
    BOOST_FOREACH( remome::remo_map::value_type v, *( _remome->_remos ) )
    {
        if( *( v.first->_centre_sequence ) == *( seq_id ) )
        {
            module::list_ptr modules = v.second;
            BOOST_FOREACH( module::ptr remo, *( modules ) )
            {
                remo_sequence::list_ptr sequences = remo->get_sequences( v.first->_centre_sequence );
                BOOST_ASSERT( sequences->size() == 1 );
                if( loc == ( *sequences )[0]->_location )
                {
                    return remo_locator( remo, v.first );
                }
            }
        }
    }

    throw std::logic_error( BIOPSY_MAKE_STRING( "Could not find remo for id: " << id ) );
}

sequence_vec_ptr
get_sequences_for_remo(
    aligned_sequence_set::ptr aligned_seqs,
    module::ptr _remo )
{
    sequence_vec_ptr result( new sequence_vec );

    result->push_back( _remo->get_sequence_for( aligned_seqs->_centre_sequence ) );
    BOOST_FOREACH( module::sequence_remo_map::value_type v, _remo->_sequences )
    {
        alignment_sequence_id::ptr seq_id = v.first;
        if( seq_id == aligned_seqs->_centre_sequence )
        {
            continue; //have already appended centre sequence
        }
        result->push_back( _remo->get_sequence_for( seq_id ) );
    }
    return result;
}



} //namespace remo
} //namespace biopsy
