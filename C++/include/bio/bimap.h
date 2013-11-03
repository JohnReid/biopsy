
#ifndef BIO_BIMAP_H_
#define BIO_BIMAP_H_

#include "bio/defs.h"

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4312)
#pragma warning(disable : 4311)
#pragma warning(pop)
#endif //_MSC_VER

BIO_NS_START


template<
	typename FromType,
	typename ToType
>
struct bimap
{
	/** tag for accessing "from" side of a bidirectional map */
	struct from{};

	/** tag for accessing "to" side of a bidirectional map */
	struct to{};

	/** The value of the things in the the map. */
	typedef std::pair< FromType, ToType > value_type;

	/* A bidirectional map can be simulated as a multi_index_container
	* of pairs of (FromType,ToType) with two indices, one
	* for each member of the pair.
	*/
	typedef boost::multi_index_container<
		value_type,
		boost::multi_index::indexed_by<
			boost::multi_index::ordered_unique<
				boost::multi_index::tag< from >,
				boost::multi_index::member<
					value_type, 
					FromType, 
					&value_type::first > >,
			boost::multi_index::ordered_unique<
				boost::multi_index::tag< to >,
				boost::multi_index::member<
					value_type,
					ToType,
					&value_type::second > >
		>
	> unique;

	/* A bidirectional multimap can be simulated as a multi_index_container
	* of pairs of (FromType,ToType) with two indices, one
	* for each member of the pair.
	*/
	typedef boost::multi_index_container<
		value_type,
		boost::multi_index::indexed_by<
			boost::multi_index::ordered_non_unique<
				boost::multi_index::tag< from >,
				boost::multi_index::member<
					value_type, 
					FromType, 
					&value_type::first > >,
			boost::multi_index::ordered_non_unique<
				boost::multi_index::tag< to >,
				boost::multi_index::member<
					value_type,
					ToType,
					&value_type::second > >
		>
	> multi;
};


BIO_NS_END

#endif //BIO_BIMAP_H_
