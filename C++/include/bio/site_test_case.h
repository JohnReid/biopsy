
#ifndef BIO_SITE_TEST_CASE_H_
#define BIO_SITE_TEST_CASE_H_

#include "bio/defs.h"
#include "bio/common.h"
#include "bio/sequence.h"
#include "bio/singleton.h"

#include <string>
#include <vector>


BIO_NS_START



/** Contains the information about a site in transfac that has been found in a remo extraction. */ 
struct SiteTestCase
{
	typedef SiteTestCase this_t;
	typedef boost::shared_ptr< this_t > ptr_t;
	typedef std::vector< ptr_t > vec_t;

	std::string name;
	TableLink site;
	seq_t centre_sequence;
	SeqList sequences; /**< The conserved sequences. */

protected:
	friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & name;
        ar & site;
        ar & centre_sequence;
        ar & sequences;
    }
};



struct SiteTestCases
	: SiteTestCase::vec_t
	, Singleton< SiteTestCases >
{
	void init_singleton();

	void init_from_biobase();

protected:
	friend class boost::serialization::access;
    template< typename Archive>
    void serialize( Archive & ar, const unsigned int version )
    {
        ar & boost::serialization::base_object< SiteTestCase::vec_t >( *this );
    }
};



BIO_NS_END

#endif //BIO_SITE_TEST_CASE_H_

