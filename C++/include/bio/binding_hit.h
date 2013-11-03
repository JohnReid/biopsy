
#ifndef BIO_BINDING_HIT
#define BIO_BINDING_HIT

#include "bio/defs.h"


#include <numeric>



BIO_NS_START



/**
One hit at a particular position in a sequence.
*/
template< typename B >
struct BindingHit
	: boost::less_than_comparable< BindingHit< B > >
	, boost::equality_comparable< BindingHit< B > >
{
#ifdef _MSC_VER
	typedef typename B binder_t;
#else //_MSC_VER
	typedef B binder_t;
#endif

	binder_t * binder;
	double p_binding;
	int position;
	unsigned length;
	bool complementary;

	BindingHit(
		binder_t * binder = 0,
		double p_binding = 0.0,
		int position = 0,
		unsigned length = 0,
		bool complementary = false )
		: binder( binder )
		, p_binding( p_binding )
		, position( position )
		, length( length )
		, complementary( complementary )
	{
	}

	const binder_t * get_binder() const
	{
		return binder;
	}

	double get_p_binding() const
	{
		return p_binding;
	}

	int get_position() const
	{
		return position;
	}

	unsigned get_length() const
	{
		return length;
	}

	void set_length( unsigned new_length )
	{
		length = new_length;
	}

	bool is_complementary() const
	{
		return complementary;
	}

	int get_end() const
	{
		return position + length;
	}

	bool operator==( const BindingHit< binder_t > & rhs ) const
	{
		return
			position == rhs.position
			&& complementary == rhs.complementary
			&& length == rhs.length
			&& boost::test_tools::check_is_close(
				p_binding,
				rhs.p_binding,
				BIO_FPC_NS::percent_tolerance( 0.01 )
			)
			&& binder == rhs.binder
			;
	}

	bool operator<( const BindingHit< binder_t > & rhs ) const
	{
		if( position < rhs.position ) return true;
		else if( ! ( rhs.position < position ) ) {
			if( complementary < rhs.complementary ) return true;
			else if( !( rhs.complementary < complementary ) ) {
				if( length < rhs.length ) return true;
				else if( ! ( rhs.length < length ) ) {
					if( p_binding < rhs.p_binding ) return true;
					else if( ! ( rhs.p_binding < p_binding ) ) {
						return binder < rhs.binder;
					}
				}
			}
		}
		return false;
	}

private:
    friend class boost::serialization::access;

	template< typename Archive >
	void serialize( Archive & ar, const unsigned int version )
	{
		ar
			& binder
			& p_binding
			& position
			& length
			& complementary;
	}
};




/**
A set of binding hits indexed by various criteria.
*/
template< typename B >
struct BindingHitSet
{
	struct position { };
	struct binder { };
	struct prob { };

	typedef boost::multi_index_container<
		BindingHit< B >,
		boost::multi_index::indexed_by<

			// ordered by position (and other members)
			boost::multi_index::ordered_non_unique<
				boost::multi_index::tag< position >,
				boost::multi_index::identity< BindingHit< B > >
			>,

			// ordered by binder
			boost::multi_index::ordered_non_unique<
				boost::multi_index::tag< binder >,
				boost::multi_index::member<
					BindingHit< B >,
					B *,
					&BindingHit< B >::binder
				>
			>,

			// ordered by binding prob
			boost::multi_index::ordered_non_unique<
				boost::multi_index::tag< prob >,
				boost::multi_index::member<
					BindingHit< B >,
					double,
					&BindingHit< B >::p_binding
				>
			>
		>
	> type;
};


template< typename HitIt >
double
get_expected_num_hits( HitIt hits_begin, HitIt hits_end )
{
	using namespace boost;

	typedef typename HitIt::value_type hit_t;

	return
		std::accumulate(
			make_transform_iterator( hits_begin, bind( &hit_t::get_p_binding, _1) ),
			make_transform_iterator( hits_end, bind( &hit_t::get_p_binding, _1) ),
			0.0 );
}


template< typename B >
double get_p_binding( const typename BindingHitSet< B >::type & hits, B * binder )
{
	using namespace boost::lambda;

	typedef B binder_t;
	typedef BindingHit< binder_t > hit_t;
	typedef typename BindingHitSet< B >::type hit_set_t;
	typedef typename hit_set_t::template index< typename BindingHitSet< binder_t >::binder >::type by_binder_t;

#ifdef _MSC_VER
	const by_binder_t & by_binder = hits.get< typename BindingHitSet< binder_t >::binder >();
#else //_MSC_VER
	const typename hit_set_t::template index< typename BindingHitSet< binder_t >::binder >::type &
		by_binder =
			::boost::multi_index::get<
				typename BindingHitSet< binder_t >::binder
			>( hits );
#endif //_MSC_VER

	double p_doesnt_bind = 1.0;
	std::pair< typename by_binder_t::const_iterator, typename by_binder_t::const_iterator > range =
		by_binder.equal_range( binder );
	for( typename by_binder_t::const_iterator mh = range.first;
		range.second != mh;
		++mh )
	{
		BOOST_ASSERT( binder == mh->get_binder() );
		p_doesnt_bind *= ( 1.0 - mh->get_p_binding() );
	}
	const double p_does_bind = 1.0 - p_doesnt_bind;

	return p_does_bind;
}


template< typename B >
std::ostream &
operator<<( std::ostream & os, const BindingHit< B > & hit )
{
	os
		<< "("
		<< hit.position
		<< ","
		<< hit.length
		<< ","
		<< ( hit.complementary ? "true" : "false" )
		<< ","
		<< hit.p_binding
		<< ","
		<< ( 0 == hit.binder ? "<no binder>" : hit.binder->get_name() )
		<< ")"
		;

	return os;
}



template< typename B, typename Archive >
void
serialise(
	Archive & ar,
	const typename BindingHitSet< B >::type & t )
{
	unsigned n = t.size();
	ar << n;
	for( typename BindingHitSet< B >::type::const_iterator h = t.begin();
		t.end() != h;
		++h )
	{
		ar << *h;
	}
}

template< typename B, typename Archive >
void
deserialise(
	Archive & ar,
	typename BindingHitSet< B >::type & t )
{
	unsigned n;
	ar >> n;
	t.clear();
	while( 0 != n )
	{
		typename BindingHitSet< B >::type::value_type value;
		ar >> value;
		t.insert( value );
		--n;
	}
}


template< typename B, bool binary >
void
serialise(
	const typename BindingHitSet< B >::type & t,
	const boost::filesystem::path & file)
{
	typedef B binder_t;

	if( binary )
	{
		boost::filesystem::ofstream stream( file, std::ios::binary );
		boost::archive::binary_oarchive ar( stream );
		serialise< binder_t >( ar, t );
	}
	else
	{
		boost::filesystem::ofstream stream( file );
		boost::archive::text_oarchive ar( stream );
		serialise< binder_t >( ar, t );
	}
}



template< typename B, bool binary >
void
deserialise(
	typename BindingHitSet< B >::type & t,
	const boost::filesystem::path & file)
{
	typedef B binder_t;

	if( binary )
	{
		boost::filesystem::ifstream stream( file, std::ios::binary );
		boost::archive::binary_iarchive ar( stream );
		deserialise< binder_t >( ar, t );
	}
	else
	{
		boost::filesystem::ifstream stream( file );
		boost::archive::text_iarchive ar( stream );
		deserialise< binder_t >( ar, t );
	}
}




BIO_NS_END



#endif //BIO_BINDING_HIT
