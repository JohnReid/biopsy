

#ifndef BIO_EQUIVALENT_FACTORS_H_
#define BIO_EQUIVALENT_FACTORS_H_

#include "bio/defs.h"
#include "bio/equivalence_partition.h"
#include "bio/singleton.h"
#include "bio/factor.h"

#include <boost/shared_ptr.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/tuple/tuple.hpp>


BIO_NS_START

struct FactorEquivalence
{
	double synonym_proportion; /**< The proportion of matching synonyms in order to have equivalence. */

	FactorEquivalence(double synonym_proportion = 0.5)
		: synonym_proportion(synonym_proportion)
	{
	}

	static Factor * get_factor(unsigned factor_acc_number);

	template <typename It>
	bool operator()(unsigned factor_acc_number, It partition_begin, It partition_end) const
	{
		Factor * factor = get_factor(factor_acc_number);
		if (0 == factor)
		{
			throw std::invalid_argument("Null factor pointer");
		}

		//is its name a synonym for one of the factors already?
		for (It p = partition_begin; partition_end != p; ++p)
		{
			if (0 == *p)
			{
				throw std::logic_error( "Null factor pointer already in partition" );
			}

			const std::string & name = factor->get_name();

			Factor * f = get_factor((*p));
			if (name == f->get_name() || f->synonyms.find(name) != f->synonyms.end())
			{
				return true;
			}
		}

		//how many of its synonyms are synonyms of the factors in the partition?
		unsigned num_synonyms_matched = 0;
		for (Factor::synonym_set_t::const_iterator s = factor->synonyms.begin();
			factor->synonyms.end() != s;
			++s)
		{
			for (It p = partition_begin; partition_end != p; ++p)
			{
				Factor * f = get_factor((*p));
				if (f->synonyms.find(*s) != f->synonyms.end())
				{
					++num_synonyms_matched;
					break;
				}
			}
		}
		if (double(num_synonyms_matched) > synonym_proportion * double(factor->synonyms.size())) //say true if over a given proportion match
		{
			return true;
		}

		return false;
	}

protected:
	friend class boost::serialization::access;
	template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
		ar & synonym_proportion;
    }
};

/**
A partition of pointers to factors into equivalence classes. The factors are keyed by their accession numbers.
*/
struct EquivalentFactors
	: EquivalencePartition< unsigned, FactorEquivalence >
	, Singleton< EquivalentFactors >
{
public:
	typedef boost::shared_ptr< EquivalentFactors > ptr_t;
	typedef EquivalencePartition< unsigned, FactorEquivalence > base_t;
	typedef boost::tuples::tuple< partition_ptr_t, partition_ptr_t > pair_t;

	/** Get the name for this partition. */
	static std::string get_name_for(partition_ptr_t);

	/** Get an accession id that identifies this factor. */
	static unsigned get_indicative_acc_id(partition_ptr_t);

	/** Create from data in biobase. */
	static ptr_t construct_from_biobase();

	/** Initialise the singleton. */
	void init_singleton();

	/** Get the factor partitions that a pssm can represent. */
	partition_set_t get_factors_for(BiobaseTablePssmEntry * pssm) const;

	/** How to initialise if cannot deserialise from disk. */
	void init_from_biobase();

	/** Get the default archive file path. */
	boost::filesystem::path get_archive_file() const;

protected:
	friend class boost::serialization::access;
	template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
		ar & boost::serialization::base_object< base_t >(*this);
    }
};



/**
Converts a variety of types to their factor partition.
*/
struct EquivalentFactorKeyTransformer
{
	typedef EquivalentFactors::partition_ptr_t result_type;

	const result_type & operator()( const result_type & partition ) const;
	result_type operator()( const TableLink & factor ) const;
};




typedef std::set<BiobaseTablePssmEntry *> pssm_set;

void get_all_pssms(pssm_set & pssms);


BIO_NS_END

BOOST_CLASS_TRACKING( BIO_NS::EquivalentFactors::partition_t, boost::serialization::track_always )

#endif //BIO_EQUIVALENT_FACTORS_H_

