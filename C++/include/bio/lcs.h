/**
@file

Copyright John Reid 2006

*/

#ifndef BIO_LCS_H_
#define BIO_LCS_H_

#ifdef _MSC_VER
# pragma once
#endif //_MSC_VER



#include <bio/defs.h>

#include <boost/multi_array.hpp>


BIO_NS_START





template<
    typename ValueT, 
    typename CharT, 
    typename CharExtractor, 
    typename StartExtractor, 
    typename EndExtractor, 
    typename ScoreExtractor >
struct LCS
{
    typedef CharT character;
    typedef std::set< character > char_set;
    typedef std::vector< character > string;

    typedef ValueT value;
    typedef std::vector< value > seq;
    typedef std::vector< seq > seq_vec;

    typedef CharExtractor char_extractor;
    typedef StartExtractor start_extractor;
    typedef EndExtractor end_extractor;
    typedef ScoreExtractor score_extractor;

    //typedef std::set< int > pos_set;
    typedef std::vector< int > pos_vec;
    typedef std::vector< pos_vec > pos_array;

    /** Index into our storage. */
    typedef std::vector< unsigned > index;

    /** Comparision operator to sort values by end position. */
    struct end_cmp
    {
        end_cmp( end_extractor end_ext ) : _end_ext( end_ext ) { }
        end_extractor _end_ext;
        bool operator()( const value & v1, const value & v2 ) const
        {
            return _end_ext( v1 ) < _end_ext( v2 );
        }
    };

    struct value_holder
    {
        typedef std::vector< value_holder > vec;
        typedef boost::shared_ptr< value_holder > ptr;

        value_holder(
            int start,
            int end,
            character c,
            double score )
            : _start( start )
            , _end( end )
            , _c( c )
            , _score( score )
        {
        }

        int _start;
        int _end;
        character _c;
        double _score;
    };

    /** Holds details of best LCS so far. */
    struct best_lcs
    {
    protected:
        best_lcs * _previous;

    public:
        typename value_holder::ptr _value;


        typedef std::vector< best_lcs > vec;

        best_lcs()
            : _previous( 0 )
        { }

        best_lcs(
            best_lcs * previous,
            double score,
            int start,
            int end,
            character c )
            : _previous( previous )
            , _value( new value_holder( start, end, c, score ) )
        { 
            BOOST_ASSERT( _previous != this );
        }

        best_lcs & operator=( const best_lcs & rhs )
        {
            BOOST_ASSERT( this != rhs._previous );

            _previous = rhs._previous;
            _value = rhs._value;

            return *this;
        }

        void set_previous( best_lcs * previous )
        {
            BOOST_ASSERT( previous != this );
            _previous = previous;
        }

        double get_score() const
        {
            BOOST_ASSERT( _previous != this );
            return ( _value ? _value->_score : 0.0 ) + ( _previous ? _previous->get_score() : 0.0 );
        }

        typename value_holder::vec & fill_values( typename value_holder::vec & values ) const
        {
            if( _previous )
            {
                _previous->fill_values( values );
            }
            if( _value )
            {
                values.push_back( *_value );
            }
            return values;
        }

        typename value_holder::vec get_values() const
        {
            typename value_holder::vec result;
            fill_values( result );
            return result;
        }

        string & fill_string( string & str ) const
        {
            if( _previous )
            {
                _previous->fill_string( str );
            }
            if( _value )
            {
                str.push_back( _value->_c );
            }
            return str;
        }

        string get_string() const
        {
            string result;
            fill_string( result );
            return result;
        }
    };

    /** Iterates across multiple sequences. */
    struct multi_iterator
    {
        multi_iterator( index begin, index end )
            : _ind( begin )
            , _begin( begin )
            , _end( end )
        {
            if( begin.size() != end.size() )
            {
                throw std::invalid_argument( "begin and end sizes are different" );
            }
        }
        index _ind;
        index _begin;
        index _end;

        unsigned get_dimensions() const
        {
            return _ind.size();
        }

        bool at_end()
        {
            for( unsigned d = 0; get_dimensions() != d; ++d )
            {
                if( _ind[ d ] != _end[ d ] )
                {
                    return false;
                }
            }
            return true;
        }

        void next()
        {
            for( unsigned d = 0; get_dimensions() != d; ++d )
            {
                _ind[ d ]++;
                if( _end[ d ] == _ind[ d ] )
                {
                    _ind[ d ] = _begin[ d ];
                }
                else
                {
                    return; //we have incremented the index
                }
            }
            //set ind to be the end
            _ind = _end;
        }
    };

    seq_vec _sequences; /**< sequences of values */
    pos_array _end_positions; /**< end positions we are interested in */
    char_extractor _char_ext; /**< gets characters from values */
    start_extractor _start_ext; /**< gets start position from values */
    end_extractor _end_ext; /**< gets end position from values */
    score_extractor _score_ext; /**< gets score from values */
    char_set _universe; /**< universe of characters in each sequence */
    typename best_lcs::vec _best; /**< storage for best LCS so far */
    index _end_position; /**< An index just past the position... */
    index _end_value; /**< An index just past the last value... */

    /**
    Gets the character associated with the value.
    */
    character get_char( const value & v ) const
    {
        return _char_ext( v );
    }


    /**
    Gets the start position associated with the value.
    */
    int get_start( const value & v ) const
    {
        return _start_ext( v );
    }


    /**
    Gets the end position associated with the value.
    */
    int get_end( const value & v ) const
    {
        return _end_ext( v );
    }


    /**
    Gets the score associated with the value.
    */
    double get_score( const value & v ) const
    {
        return _score_ext( v );
    }

    unsigned get_num_seqs() const
    {
        return _sequences.size();
    }

    unsigned get_storage_size() const
    {
        unsigned storage_size = 1;
        BOOST_FOREACH( const pos_vec & positions, _end_positions )
        {
            storage_size *= positions.size();
        }
        return storage_size;
    }

    multi_iterator get_end_iterator() const
    {
        return multi_iterator( index( get_num_seqs(), 0 ), _end_position );
    }

    const value & get_value( unsigned seq_idx, unsigned char_idx ) const
    {
        return _sequences[ seq_idx ][ char_idx ];
    }

    /** Get a value iterator s.t. we range over all values with the
    specified ends... */
    multi_iterator get_value_range_for_end( const index & end_ind )
    {
        //for each dimension
        index begin( get_num_seqs() );
        index end( get_num_seqs() );
        //std::cout << "Looking for: ";
        for( unsigned d = 0; get_num_seqs() != d; ++d )
        {
            const int end_looking_for = _end_positions[ d ][ end_ind[ d ] ];
            //std::cout << end_looking_for << " ";
            bool found_first = false;
            end[ d ] = _sequences[ d ].size();
            for( unsigned i = 0; _sequences[ d ].size() != i; ++i )
            {
                const int end_here = get_end( get_value( d, i ) );
                if( ! found_first )
                {
                    if( end_here == end_looking_for )
                    {
                        found_first = true;
                        begin[ d ] = i;
                    }
                }

                if( end_here > end_looking_for )
                {
                    end[ d ] = i;
                    break;
                }
            }
        }
        //std::cout << "\n";
        return multi_iterator( begin, end );
    }


    /** Given the index get the offset in the best storage vector. */
    unsigned get_best_index( const index & ind )
    {
        unsigned result = 0;
        unsigned multiplier = 1;
        for( unsigned i = 0; _end_position.size() != i; ++i )
        {
            result += multiplier * ind[ i ];
            multiplier *= _end_position[ i ];
        }
        return result;
    }

    best_lcs & get_best( const index & ind )
    {
        return _best[ get_best_index( ind ) ];
    }

    best_lcs get_best()
    {
        if( _best.empty() )
        {
            return best_lcs();
        }

        return *( _best.rbegin() );
    }

    /**
    Get all characters that are in every sequence.
    Returns a set of those characters.
    */
    template< typename InputSeqRange >
    char_set get_universe( InputSeqRange sequences )
    {
        typedef typename boost::range_value< InputSeqRange >::type input_seq;
        using namespace boost;

        char_set result;

        //for each sequence
        bool is_first_seq = true;
        BOOST_FOREACH( const input_seq & _input_seq, sequences )
        {
            //get the set of all characters...
            char_set set_for_seq;
            BOOST_FOREACH( const value & v, _input_seq )
            {
                set_for_seq.insert( get_char( v ) );
            }

            //if this is the first sequence take it as our universe so far
            if( is_first_seq )
            {
                result = set_for_seq;
                is_first_seq = false;
            }
            else
            {
                char_set tmp = result;
                result.clear();
                std::set_intersection(
                    set_for_seq.begin(),
                    set_for_seq.end(),
                    tmp.begin(),
                    tmp.end(),
                    std::inserter( result, result.begin() ) );
            }
        }

        return result;
    }


    /**
    Build the sequences of values. Ignore those values whose characters are not in
    the character universe set.
    */
    template< typename InputSeqRange >
    void build_sequences( InputSeqRange sequences )
    {
        //build the sequence vectors from those values whose chars are in the universe
        //also build a vector of sets of end positions
        typedef typename boost::range_value< InputSeqRange >::type input_seq;
        BOOST_FOREACH( const input_seq & _input_seq, sequences )
        {
            seq _seq;

            BOOST_FOREACH( const value & v, _input_seq )
            {
                if( _universe.end() != _universe.find( get_char( v ) ) )
                {
                    _seq.push_back( v );
                }
            }

            //sort the sequence by the end position of each value
            std::sort( _seq.begin(), _seq.end(), end_cmp( _end_ext ) );

            _sequences.push_back( _seq );
            _end_value.push_back( _seq.size() );
        }
    }

    /**
    For each sequence calculate the set of all end values.
    */
    void build_end_values()
    {
        BOOST_FOREACH( const seq & _s, _sequences )
        {
            pos_vec end_positions;

            BOOST_FOREACH( const value & v, _s )
            {
                int end = get_end( v );
                if( end_positions.end() == std::find( end_positions.begin(), end_positions.end(), end ) )
                {
                    end_positions.push_back( end );
                }
            }

            _end_positions.push_back( end_positions );
        }

        BOOST_FOREACH( const pos_vec & positions, _end_positions )
        {
            _end_position.push_back( positions.size() );
        }
    }

    /**
    Resize the array that holds the best LCS's up to each combination of end positions.
    */
    void
    allocate_storage()
    {
        _best.resize( get_storage_size() );
    }


    /**
    Look at all the combinations of values that have the end positions indexed by ind.
    If a combination has the same character in each value then check to see if we have
    superceded the previous best LCS. If so update...
    */
    void update_best_if_common( const index & ind, best_lcs & previous_best )
    {
        if( ind.empty() )
        {
            return;
        }

        //do we have a common character
        character c = get_char( get_value( 0, ind[ 0 ] ) );
        for( unsigned d = 1; ind.size() != d; ++d )
        {
            if( get_char( get_value( d, ind[ d ] ) ) != c )
            {
                //not common character
                return;
            }
        }

        //find out the best lcs up to the point where the common character starts

        //where does the common character start and what is the average score?
        pos_vec char_starts;
        double s = 0.0;
        for( unsigned d = 0; ind.size() != d; ++d )
        {
            const value & v = get_value( d, ind[ d ] );
            s += get_score( v );
            char_starts.push_back( get_start( v ) );
        }
        s /= ind.size();

        //find the last end point before or at the start
        index previous_ind;
        bool found_previous_end_point = true;
        for( unsigned d = 0; ind.size() != d; ++d )
        {
            unsigned i = 0;
            for( ; _end_positions[ d ].size() != i; ++i )
            {
                if( _end_positions[ d ][ i ] > char_starts[ d ] )
                {
                    break;
                }
            }
            //was the start before the first end point?
            if( 0 == i )
            {
                found_previous_end_point = false;
                break;
            }
            previous_ind.push_back( i - 1 );
        }

        best_lcs * lhs = 0;
        double lhs_score = 0.0;
        //did we find a left hand side?
        if( found_previous_end_point )
        {
            lhs = boost::addressof( get_best( previous_ind ) );
            lhs_score = lhs->get_score();
        }
        //have we done better than our previous best?
        if( s + lhs_score > previous_best.get_score() )
        {
            const value & v = get_value( 0, ind[ 0 ] );
            BOOST_ASSERT( lhs != boost::addressof( previous_best ) );
            previous_best =
                best_lcs(
                    lhs,
                    s,
                    get_start( v ),
                    get_end( v ),
                    c );
        }
    }


    /**
    Solve the sub-problem of calculating the LCS up to the end position indexed by end_it.
    */
    void
    calculate_best( multi_iterator & end_it )
    {
        best_lcs b;
        best_lcs & to_assign_to = get_best( end_it._ind );

        //look one step back in each dimension to see which is best
        for( unsigned d = 0; end_it._ind.size() != d; ++d )
        {
            if( 0 == end_it._ind[ d ] )
            {
                //at beginning of this dimension
                continue;
            }

            //check if one before is better...
            end_it._ind[ d ]--;
            best_lcs & prev_best = get_best( end_it._ind );
            const double prev_best_score = prev_best.get_score();
            const double score = b.get_score();
            if( prev_best_score > score )
            {
                b.set_previous( boost::addressof( prev_best ) );
            }
            end_it._ind[ d ]++;
        }

        //now check to see if we can add a common character to improve on carried forward best subsequences
        //so we need to iterate through all combinations of values from all the sequences where the end
        //position is the end position we are looking at...
        for( multi_iterator value_it = get_value_range_for_end( end_it._ind ); ! value_it.at_end(); value_it.next() )
        {
            update_best_if_common( value_it._ind, b );
        }

        //assign best to storage
        to_assign_to = b;
    }


    /**
    Calculate the best LCS over all combinations of end points.
    */
    void
    calculate_best()
    {
        //make sure we have enough room to store calculations
        allocate_storage( );

        //for each end index, solve the sub-problem
        for( multi_iterator end_it = get_end_iterator(); ! end_it.at_end(); end_it.next() )
        {
            calculate_best( end_it );
        }
    }

    template< typename InputSeqRange >
    LCS(
        InputSeqRange sequences,
        char_extractor char_ext = char_extractor(),
        start_extractor start_ext = start_extractor(),
        end_extractor end_ext = end_extractor(),
        score_extractor score_ext = score_extractor() )
        : _char_ext( char_ext )
        , _start_ext( start_ext )
        , _end_ext( end_ext )
        , _score_ext( score_ext )
    {
        _universe = get_universe( sequences );

        build_sequences( sequences );

        build_end_values( );

    };
};






BIO_NS_END

#endif //BIO_LCS_H_

