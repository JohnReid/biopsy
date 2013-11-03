
#ifndef BIO_BIOBASE_PROXY_H_
#define BIO_BIOBASE_PROXY_H_

#include "bio/defs.h"

BIO_NS_START




template< typename Retriever >
struct Proxy
{
	typedef Retriever retriever_t;
	typedef typename retriever_t::argument_type index_t;
	typedef typename retriever_t::result_type object_t;

	retriever_t retriever;
	index_t index;

	Proxy( retriever_t retriever, index_t index )
		: retriever( retriever )
		, index( index )
	{
	}

	object_t operator()() const
	{
		return retriever( index );
	}
};




BIO_NS_END

#define BIO_BIOBASE_PROXY_H_
