
#ifndef BIO_RAII_H_
#define BIO_RAII_H_

#include "bio/defs.h"

BIO_NS_START

/** Implements "resource acquisition is initialisation" idiom. Users need to specialise destructor to free resource. */
template <class Resource>
struct RAII
{
	Resource resource;

	RAII(Resource resource);

	/** Needs to be specialised to free the resource. */
	~RAII();

	/** Const cast operator. */
	operator const Resource () const;

	/** Cast operator. */
	operator Resource ();
};

template <class Resource>
RAII<Resource>::RAII(Resource resource)
: resource(resource)
{
}

template <class Resource>
RAII<Resource>::~RAII()
{
    // static assert seems to always trigger now. Use a runtime exception instead.
    // BOOST_STATIC_ASSERT("Specialise this to free your resource.");
    throw std::logic_error( "Specialise this to free your resource." );
}

template <class Resource>
RAII<Resource>::operator const Resource() const
{
	return resource;
}

/** Cast operator. */
template <class Resource>
RAII<Resource>::operator Resource ()
{
	return resource;
}


BIO_NS_END


#endif //BIO_RAII_H_
