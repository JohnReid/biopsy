#ifndef BIO_SCORE_MAP_H_
#define BIO_SCORE_MAP_H_

#include "bio/defs.h"

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4312)
#pragma warning(disable : 4311)
#endif //_MSC_VER

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/mem_fun.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_io.hpp>

#ifdef _MSC_VER
#pragma warning(pop)
#endif //_MSC_VER

#include <iostream>


BIO_NS_START




/** Maps elements to scores and vice versa. */
template< typename T >
struct ScoreMap
{
public:
	struct key{};
	struct score{};

public:
	struct value_type
	{
		T key;
		double score;

		value_type(const T & key, double score)
			: key(key), score(score)
		{
		}
	};

	typedef boost::multi_index::template multi_index_container<
		value_type,
		boost::multi_index::template indexed_by<
			boost::multi_index::ordered_unique<
				boost::multi_index::template tag< key >,
				boost::multi_index::template member< value_type, T, &value_type::key >
			>,
			boost::multi_index::ordered_non_unique<
				boost::multi_index::template tag< score >,
				boost::multi_index::template member< value_type, double, &value_type::score >
			>
		>
	> type;

	typedef typename type::template index< key >::type key_index;
	typedef typename type::template index< score >::type score_index;
};






BIO_NS_END

#endif //BIO_SCORE_MAP_H_
