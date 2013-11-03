#ifndef BIO_SERIALISABLE_H_
#define BIO_SERIALISABLE_H_

#include "bio/defs.h"
#include "bio/environment.h"

#include <boost/filesystem/path.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/concept_check.hpp>
#include <boost/bind.hpp>
#include <boost/timer.hpp>



BIO_NS_START




template< bool binary, typename T>
void
serialise(
	const T & t,
	const boost::filesystem::path & file)
{
	if( binary )
	{
		boost::filesystem::ofstream stream(file, std::ios::binary);
		boost::archive::binary_oarchive(stream) << t;
	}
	else
	{
		boost::filesystem::ofstream stream(file);
		boost::archive::text_oarchive(stream) << t;
	}
}



template< bool binary, typename T >
void
deserialise(
	T & t,
	const boost::filesystem::path & file)
{
	if( binary )
	{
		boost::filesystem::ifstream stream(file, std::ios::binary);
		boost::archive::binary_iarchive(stream) >> t;
	}
	else
	{
		boost::filesystem::ifstream stream(file);
		boost::archive::text_iarchive(stream) >> t;
	}
}



template< bool binary, typename T >
boost::shared_ptr< T >
deserialise(
	const boost::filesystem::path & file)
{
	boost::shared_ptr< T > result(new T);

	deserialise< binary, T >( *result, file );

	return result;
}



/**
Try to deserialise from disk. If not print error and return false.
*/
template<
	bool binary,				// binary serialisation format?
	typename T					// the type of object we are serialising
>
bool try_to_deserialise(
	T & object,
	const boost::filesystem::path & archive_file )
{
	typedef T object_t;

	//try to deserialise
	bool deserialised = false;
	try
	{
		if( ! boost::filesystem::exists( archive_file ) )
		{
			*( BioEnvironment::singleton().get_log_stream() ) << "\"" << archive_file._BOOST_FS_NATIVE() << "\" does not exist\n";
		}
		else
		{
			boost::timer timer;
			BIO_NS::deserialise< binary, object_t >( object, archive_file );
			*( BioEnvironment::singleton().get_log_stream() ) << "Deserialised \"" << archive_file._BOOST_FS_NATIVE() << "\" - " << timer.elapsed() << "s\n";

			deserialised = true;
		}
	}
	catch( const std::exception & exception )
	{
		*( BioEnvironment::singleton().get_log_stream() )
			<< exception.what()
			<< ": Could not deserialise \""
			<< archive_file._BOOST_FS_NATIVE()
			<< "\"\n";
	}
	catch(...)
	{
		*( BioEnvironment::singleton().get_log_stream() )
			<< "Could not deserialise \""
			<< archive_file._BOOST_FS_NATIVE()
			<< "\"\n";
	}

	return deserialised;
}


/**
Deserialise from disk. If not possible init using Init functor.
*/
template<
	bool binary,				// binary serialisation format?
	typename T,					// the type of object we are serialising
	typename Init				// function object to initialise object if deserialising fails
>
void deserialise_or_init(
	T & object,
	const boost::filesystem::path & archive_file,
	Init init )
{
	typedef T object_t;

	if( ! try_to_deserialise< binary >( object, archive_file ) )
	{
		boost::timer timer;

		//if we couldn't deserialise, construct from scratch
		init( object );
		*( BioEnvironment::singleton().get_log_stream() ) << "Initialised \"" << archive_file._BOOST_FS_NATIVE() << "\" - " << timer.elapsed() << "s\n";
		timer.restart();

		//save for next time
		serialise< binary, object_t >( object, archive_file );
		*( BioEnvironment::singleton().get_log_stream() ) << "Serialised \"" << archive_file._BOOST_FS_NATIVE() << "\" - " << timer.elapsed() << "s\n";
	}
}






BIO_NS_END

#endif //BIO_SERIALISABLE_H_

