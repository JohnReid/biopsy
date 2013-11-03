/* Copyright John Reid 2007
*/

#include "bio-pch.h"


#include "bio/defs.h"



#include "bio/random.h"

#include <boost/date_time/posix_time/posix_time_types.hpp>
#include <boost/date_time/posix_time/time_parsers.hpp>
using namespace boost;
using namespace boost::posix_time;

#include <limits>
#include <iostream>
using namespace std;


#undef min
#undef max


BIO_NS_START

default_rng_t default_rng;

namespace impl
{
	typedef boost::uniform_int<size_t> dist_t;
	typedef boost::variate_generator<default_rng_t, dist_t> var_gen_t;

	dist_t dist(numeric_limits<size_t>::min(), numeric_limits<size_t>::max());
	var_gen_t var_gen(default_rng, dist);

	bool default_rng_already_seeded = false;
};

double
get_uniform_01()
{
	ensure_default_rng_seeded();

	typedef boost::uniform_real<> distribution_type;
	typedef boost::variate_generator<default_rng_t &, distribution_type> gen_type;
	return gen_type(default_rng, distribution_type(0,1))();
}

size_t
get_uniform_index(size_t max)
{
	if (0 == max)
	{
		throw std::logic_error( "max must be >= 1" );
	}

	// Define a uniform random number distribution of integer values between
	// 0 and max-1 inclusive.
	typedef boost::uniform_int<size_t> distribution_type;
	typedef boost::variate_generator<default_rng_t &, distribution_type> gen_type;

	return gen_type(default_rng, distribution_type(0, max - 1))();
}


size_t get_random_seed()
{
#ifdef _DEBUG
	return 1; //guarantee reproducibility
#else
	const ptime now(posix_time::second_clock::universal_time());
	const ptime then(posix_time::from_iso_string("19720112T000000"));
	const time_duration diff = now - then;
	return (size_t) diff.total_seconds();
#endif
}

void
seed_default_rng(size_t seed)
{
	cout << "Using random seed: " << seed << endl;

	//seed the random number generator
	default_rng.seed((default_rng_t::result_type) seed);

	impl::default_rng_already_seeded = true;
}

void
ensure_default_rng_seeded()
{
	if (! impl::default_rng_already_seeded)
	{
		seed_default_rng(get_random_seed());
	}
}

BIO_NS_END
