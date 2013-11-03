

#ifdef _MSC_VER
#pragma warning( push )
#pragma warning( disable : 4099 )
#endif // _MSC_VER
# include <boost/archive/text_oarchive.hpp>
# include <boost/archive/text_iarchive.hpp>
# include <boost/serialization/vector.hpp>
# include <boost/serialization/map.hpp>
# include <boost/serialization/export.hpp>
# include <boost/shared_ptr.hpp>
# include <boost/test/unit_test.hpp>
# include <boost/assign/list_of.hpp>
#ifdef _MSC_VER
#pragma warning( pop )
#endif // _MSC_VER


#include <string>
#include <iostream>
#include <sstream>
#include <map>
using namespace std;


#include <boost/test/utils/wrap_stringstream.hpp>
#define MAKE_STRING(x) (boost::wrap_stringstream().ref() << x).str()


//#define VERBOSE_CHECKING




namespace my_ns
{

struct AbstractModel
{
	typedef boost::shared_ptr< AbstractModel > ptr_t;

	struct parameter_t
	{
		virtual ~parameter_t() { }
		virtual AbstractModel * get_model() const = 0;
	};

	virtual ~AbstractModel() { }
	virtual std::string get_name() const = 0;
	virtual const parameter_t * get_parameters() const = 0;
};


template< typename Key >
struct DerivedModel
	: AbstractModel
{
	struct parameter_t
		: AbstractModel::parameter_t
	{
		Key key;

		parameter_t( const Key & key = Key() )
			: key( key )
		{
		}

		virtual AbstractModel * get_model() const
		{
			return new DerivedModel( key ); //in my code this is done by a singleton cache object...
		}

	private:
		friend class boost::serialization::access;
		template< typename Archive >
		void serialize( Archive & ar, const unsigned int version )
		{
			boost::serialization::void_cast_register<
				typename DerivedModel::parameter_t,
				typename AbstractModel::parameter_t >( 0, 0 );

			ar & key;
		}
	};

	parameter_t param;

	DerivedModel( const Key & key = Key() )
		: param( key )
	{
	}

	virtual std::string get_name() const
	{
		return MAKE_STRING( "DerivedModel: " << param.key );
	}

	virtual const AbstractModel::parameter_t * get_parameters() const
	{
		return &param;
	}
};





struct Datum
{
	AbstractModel * value;
	const AbstractModel * const_value;

	Datum( AbstractModel * datum = 0 )
		: value( datum )
		, const_value( datum )
	{
	}

	bool operator==( const Datum & rhs ) const
	{
		return value->get_name() == rhs.value->get_name();
	}

private:
    friend class boost::serialization::access;
    template< typename Archive >
    void serialize( Archive & ar, const unsigned int version )
    {
		ar & value;
		ar & const_value;
    }
};

std::ostream &
operator<<( std::ostream & os, const Datum & datum )
{
	return os << datum.value->get_name();
}

} //my_ns




template< typename is_saving >
struct do_abstract_model_serialisation
{
};

//partially specialise for saving archive
template<  >
struct do_abstract_model_serialisation< boost::mpl::bool_< true > >
{
	template< typename Archive >
	void operator()( Archive & ar, const my_ns::AbstractModel * & abs ) const
	{
		const my_ns::AbstractModel::parameter_t * const parameters = abs->get_parameters();
		ar << parameters;
	}

	template< typename Archive >
	void operator()( Archive & ar, my_ns::AbstractModel * & abs ) const
	{
		const my_ns::AbstractModel::parameter_t * const parameters = abs->get_parameters();
		ar << parameters;
	}
};

//partially specialise for loading archive
template<  >
struct do_abstract_model_serialisation< boost::mpl::bool_< false > >
{
	template< typename Archive >
	void operator()( Archive & ar, my_ns::AbstractModel * & abs ) const
	{
		my_ns::AbstractModel::parameter_t * parameters;
		ar >> parameters;
		abs = parameters->get_model();
	}

	template< typename Archive >
	void operator()( Archive & ar, const my_ns::AbstractModel * & abs ) const
	{
		my_ns::AbstractModel::parameter_t * parameters;
		ar >> parameters;
		abs = parameters->get_model();
	}
};


template< typename Archive >
Archive &
operator&( Archive & ar, const typename my_ns::AbstractModel * & abs )
{
	do_abstract_model_serialisation< typename Archive::is_saving >()( ar, abs );

	return ar;
}



template< typename Archive >
Archive &
operator&( Archive & ar, typename my_ns::AbstractModel * & abs )
{
	do_abstract_model_serialisation< typename Archive::is_saving >()( ar, abs );

	return ar;
}

namespace boost { namespace archive {
template< typename Archive >
inline void save( Archive & ar, const typename my_ns::AbstractModel * &t )
{
}
} } //archive, boost

BOOST_IS_ABSTRACT( my_ns::AbstractModel::parameter_t )
BOOST_CLASS_EXPORT( my_ns::DerivedModel< int >::parameter_t )
BOOST_CLASS_EXPORT( my_ns::DerivedModel< std::string >::parameter_t )



void
check_serialisation_strategy()
{
	cout << "******* check_serialisation_strategy()\n";

	using namespace my_ns;

	typedef std::vector< Datum > data_vec_t;
	const data_vec_t data = boost::assign::list_of
		( Datum( DerivedModel< int >::parameter_t( 1 ).get_model() ) )
		( Datum( DerivedModel< std::string >::parameter_t( "hello" ).get_model() ) )
		( Datum( DerivedModel< int >::parameter_t( 2 ).get_model() ) )
		( Datum( DerivedModel< std::string >::parameter_t( "goodbye" ).get_model() ) )
		( Datum( DerivedModel< int >::parameter_t( 1 ).get_model() ) )
		;

	const Datum datum( 0 );
	const Datum * const datum_ptr( 0 );
	{
		std::stringstream ss;
		boost::archive::text_oarchive( ss ) << datum_ptr;

		Datum * copy_of_datum;
		boost::archive::text_iarchive( ss ) >> copy_of_datum;
	}

#ifdef VERBOSE_CHECKING
	std::copy(
		data.begin(),
		data.end(),
		std::ostream_iterator< Datum >( std::cout, "\n" ) );
#endif

#if 0
	//serialise and deserialise
	{
		std::stringstream ss;
		boost::archive::text_oarchive( ss ) << data;

		data_vec_t copy_of_data;
		boost::archive::text_iarchive( ss ) >> copy_of_data;

		BOOST_REQUIRE( copy_of_data.size() == data.size() );
		BOOST_CHECK( std::equal( data.begin(), data.end(), copy_of_data.begin() ) );
	}
#endif

	typedef std::map< AbstractModel *, int > model_int_map_t;
	const model_int_map_t model_int_map = boost::assign::map_list_of
		( data[0].value, 0 )
		( data[1].value, 1 )
		( data[2].value, 2 )
		;

	//serialise and deserialise
	{
		std::stringstream ss;
		boost::archive::text_oarchive( ss ) << model_int_map;

		model_int_map_t copy_of_map;
		boost::archive::text_iarchive( ss ) >> copy_of_map;

		BOOST_REQUIRE( copy_of_map.size() == model_int_map.size() );
		BOOST_CHECK( std::equal( model_int_map.begin(), model_int_map.end(), copy_of_map.begin() ) );
	}
}



void
register_serialisation_strategy_tests( boost::unit_test::test_suite * test )
{
	test->add( BOOST_TEST_CASE( &check_serialisation_strategy ), 0);
}
