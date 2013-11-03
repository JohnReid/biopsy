

#ifndef BIO_EQUIVALENCE_PARTITION_H_
#define BIO_EQUIVALENCE_PARTITION_H_

#include "bio/defs.h"

#include <boost/shared_ptr.hpp>

#include <set>


BIO_NS_START


/** Defines a partition of objects of type C into equivalence classes, defined by an equivalence operator, Equiv. */
template <typename C, typename Equiv>
class EquivalencePartition
{
public:
	typedef C class_t; /**< The type of object we are partitioning. */
	typedef Equiv equiv_t; /**< The relationship that partitions the objects. */
	typedef std::set< class_t > partition_t; /**< One partition. */
	typedef typename partition_t::const_iterator iterator; /**< Iterator over one partition. */
	typedef boost::shared_ptr< partition_t > partition_ptr_t; /**< A pointer to one partition. */
	typedef std::set< partition_ptr_t > partition_set_t; /**< A set of partitions. */
	typedef typename partition_set_t::const_iterator partition_set_iterator; /**< Iterator over all partitions. */

protected:
	equiv_t equiv;  /**< The relationship that partitions the objects. */
	partition_set_t partitions; /** The set of partitions. */

public:
	EquivalencePartition(equiv_t equiv = equiv_t())
		: equiv(equiv)
	{
	}

	unsigned num_partitions() const
	{
		return partitions.size();
	}

	partition_set_iterator begin() const
	{
		return partitions.begin();
	}

	partition_set_iterator end() const
	{
		return partitions.end();
	}

	/** Get the partition the object is a member of. Throws exception if none. */
	partition_ptr_t get_partition(class_t c) const
	{
		partition_ptr_t result = find_partition(c);
		if (0 == result)
		{
			throw std::logic_error( "Could not find partition" );
		}

		return result;
	}

	/** Find the partition the object is a member of. Returns 0 for none. */
	partition_ptr_t find_partition(class_t c) const
	{
		for (typename partition_set_t::const_iterator p = partitions.begin();
			partitions.end() != p;
			++p)
		{
			//is the object is already in the partition
			if (p->get()->find(c) != p->get()->end())
			{
				//return the set (partition)
				return *p;
			}
		}

		//return a null pointer
		return partition_ptr_t();
	}

	/** Find the partition(s) the object would be a member of. */
	partition_set_t which_partitions(class_t c) const
	{
		partition_set_t result;

		for (typename partition_set_t::const_iterator p = partitions.begin();
			partitions.end() != p;
			++p)
		{
			//if the object is already in the partition or it would be under the equivalence relationship
			if (p->get()->find(c) != p->get()->end() || equiv(c, p->get()->begin(), p->get()->end()))
			{
				//return the set
				result.insert(*p);
			}
		}

		return result;
	}

	/** Add an object to an existing partition or create a new one. This can invalidate existing pointers to partitions. */
	partition_ptr_t add_object(class_t c)
	{
		partition_set_t equivalent_partitions = which_partitions(c);

		//did we find any partitions?
		if (0 == equivalent_partitions.size())
		{
			//no - so create new one and add
			partition_ptr_t partition(new partition_t());
			partition->insert(c);

			partitions.insert(partition);

			return partition;
		}
		else
		{
			//yes - so merge the partitions
			typename partition_set_t::iterator merge_into_partition = equivalent_partitions.begin();
			typename partition_set_t::iterator p = merge_into_partition;
			++p;
			for ( ;
				equivalent_partitions.end() != p;
				++p)
			{
				//merge the partition
				(*merge_into_partition)->insert((*p)->begin(), (*p)->end());

				//remove from our list of partitions
				partitions.erase(*p);
			}

			//add our object
			(*merge_into_partition)->insert(c);

			return *merge_into_partition;
		}
	}

protected:
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & ar, const unsigned int version)
	{
		ar & equiv;
		ar & partitions;
	}
};



BIO_NS_END


#endif //BIO_EQUIVALENCE_PARTITION_H_

