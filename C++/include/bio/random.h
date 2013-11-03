#ifndef BIO_RANDOM_H_
#define BIO_RANDOM_H_

#include "bio/defs.h"

#include <boost/random.hpp>

#include <algorithm>
#include <vector>


BIO_NS_START


typedef boost::mt19937 default_rng_t;
extern default_rng_t default_rng;


//gets a random number to seed rng's with
size_t get_random_seed();

void ensure_default_rng_seeded();
void seed_default_rng(size_t seed);

double get_uniform_01();
size_t get_uniform_index(size_t max); //returns number in [0,max-1]

struct CreateRandomUniform01Partition
{
	template <class InsIt>
	void operator()(size_t num_partitions, InsIt ins_it) const
	{
		if (0 == num_partitions) {
			throw std::logic_error( "Need positive num_partitions" );
		}

		//put random numbers in a partition
		std::vector<double> partitions;
		for (size_t i = 0; i < num_partitions - 1; ++i)
		{
			partitions.push_back(get_uniform_01());
		}

		//sort the numbers
		std::sort(partitions.begin(), partitions.end());

		//copy to the insert iterator
		std::copy(partitions.begin(), partitions.end(), ins_it);
	}
};


struct GetUniform01PartitionSizes
{
	template <class PartIt, class InsIt>
	void operator()(PartIt begin, PartIt end, InsIt ins_it) const
	{
		typedef typename PartIt::value_type value_type;

		value_type last_boundary = 0.0;
		for ( ; begin != end; ++begin) {
			*ins_it++ = *begin - last_boundary;
			last_boundary = *begin;
		}
		*ins_it++ = 1.0 - last_boundary;
	}
};


/** Generates num_partitions uniform random numbers that sum to 1. */
template <class InsIt>
void
gen_prob_dist(size_t num_partitions, InsIt ins_it)
{
	std::vector<double> partitions;

	//create the partition
	CreateRandomUniform01Partition()(num_partitions, inserter(partitions, partitions.begin()));

	//get the sizes of the partition
	GetUniform01PartitionSizes()(partitions.begin(), partitions.end(), ins_it);
}


BIO_NS_END



#endif //BIO_RANDOM_H_
