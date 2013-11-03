#ifndef BIO_PSSM_LIKELIHOOD_CACHE_H_
#define BIO_PSSM_LIKELIHOOD_CACHE_H_


#include "bio/defs.h"
#include "bio/biobase_likelihoods.h"
#include "bio/biobase_db.h"
#include "bio/singleton.h"

#include <string>



BIO_NS_START

/** Caches the likelihoods of biobase scores (given that they bind) for pssms. */
class PssmLikelihoodCache
	: public Singleton< PssmLikelihoodCache >
{
	friend struct Singleton< PssmLikelihoodCache >;

private:
    friend class boost::serialization::access;
    // When the class Archive corresponds to an output archive, the
    // & operator is defined similar to <<.  Likewise, when the class Archive
    // is a type of input archive the & operator is defined similar to >>.
    template< typename Archive >
    void serialize(Archive & ar, const unsigned int version)
	{
        ar & likelihoods;
    }

	void init_singleton();

public:
	typedef TableLink key_t;

protected:
	typedef std::map< key_t, BiobaseLikelihoods > likelihood_map_t;

	BiobaseDb & biobase_db;
	likelihood_map_t likelihoods;

public:
	/** Constructed with which biobase db we use to generate the likelihoods from. */
	PssmLikelihoodCache( BiobaseDb & biobase_db = BiobaseDb::singleton() );

	/** Gets the likelihoods for the given pssm. These will be calculated as needed. */
	const BiobaseLikelihoods * get_likelihoods( const key_t & key );

	bool operator==(const PssmLikelihoodCache & rhs) const;

	void populate_from_biobase();
};

BIO_NS_END


#endif //BIO_PSSM_LIKELIHOOD_CACHE_H_
