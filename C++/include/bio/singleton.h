#ifndef BIO_SINGLETON_H_
#define BIO_SINGLETON_H_

#include "bio/defs.h"
#include "bio/useradmin.h"

#include <boost/utility.hpp>

BIO_NS_START

/**
 * Note that ## means the NAME parameter is concatenated with the following text
 */
#define ADD_STATIC_SINGLETON_VARIABLE(TYPE,NAME) TYPE NAME;  \
	static TYPE NAME##Get() {return singleton().NAME;}; \
	static void NAME##Set(TYPE newVal){singleton().NAME = newVal;};

/**	For use when this is declared as an external Python Interface */
#define ADD_STATIC_PROPERTY(CLASS,PROP) .add_static_property( #PROP, &CLASS::PROP##Get,&CLASS::PROP##Set)

/**
 * Inherit from this class if you want a common singleton for all users.
 */
template< typename T >
struct Singleton
: boost::noncopyable
{
	static T & singleton()
	{
		typedef T object_t;
		typedef boost::shared_ptr< object_t > ptr_t;

		static ptr_t _singleton;

		if( 0 == _singleton )
		{
			_singleton.reset( new object_t() );
		
			try
			{
				_singleton->init_singleton();
			}
			catch( ... )
			{
				_singleton.reset();
				throw;
			}
		}
		return *_singleton;
	}
private:
	void init_singleton() {}

};




/**
 * This class provides for a individual 'singleton' for each user, allowing multiple
 * independent users on a multi user server.   They are identified by the username
 * that has to be passed as part of the soap message and is used to set
 * bio:user_admin::user
 */

template< typename T >
struct UserSingleton
: boost::noncopyable
{
	typedef T object_t;
	typedef boost::shared_ptr< object_t > ptr_t;

	static T & singleton()
	{

		typedef std::map< std::string, ptr_t > per_user_singleton_map;
		static per_user_singleton_map _per_user_singleton;

		const std::string username = user_admin::userGet();

		// no need for find, insert only inserts when not there already
		//typename ptrMap_t::iterator it = _singletonMap.find( username );
		std::pair< typename per_user_singleton_map::iterator, bool > insert_result =
				_per_user_singleton.insert(
				typename per_user_singleton_map::value_type(
					username,
					ptr_t( new object_t() )
				)
			);

		// if we inserted the object, we need to initialise it
		if( insert_result.second ) { // initialise if we inserted
			insert_result.first->second->init_singleton();
		}

		// return the object
		return *( insert_result.first->second );
	}

private:
	void init_singleton() {}

	/* Can't use boost::noncopyable as we need to copy the class when making
	user specific variants so declare copy constructors privately for limited access */
	/* Why do you need to do this? I have made it boost::noncopyable again. I had
	 a problem where my parameter settings were not being used properly. I think
	this UserSingleton code needs looking at. John. */
//	UserSingleton( const UserSingleton& ){};
	// why declare a do-nothing operator=?
	//const UserSingleton& operator=( const UserSingleton& ){};

};


BIO_NS_END

#endif //BIO_SINGLETON_H_
