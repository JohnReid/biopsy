

#ifndef BIO_ITERATOR_H_
#define BIO_ITERATOR_H_

#include "bio/defs.h"

#include <boost/iterator/iterator_facade.hpp>


BIO_NS_START

/** An iterator that always points at a single value. */
template < typename T >
struct single_value_iterator
  : boost::iterator_facade< single_value_iterator< T >, T, boost::forward_traversal_tag >
{
private:
	T & value;

public:
	single_value_iterator(T & value)
	  : value(value)
	{
	}

private:
	friend class boost::iterator_core_access;

	void increment()
	{
	}

	bool equal(single_value_iterator< T > const& other) const
	{
		return &(this->value) == &(other.value);
	}

	T & dereference() const { return value; }
};




BIO_NS_END

#endif //BIO_ITERATOR_H_

