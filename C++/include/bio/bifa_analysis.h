#ifndef BIO_BIFA_ANALYSIS_H_
#define BIO_BIFA_ANALYSIS_H_

#include "bio/defs.h"
#include "bio/binding_model.h"
#include "bio/sequence.h"
#include "bio/run_match.h"

#include <map>


BIO_NS_START




typedef BindingModel::hit_set_t bifa_hits_t;




struct BiFaAnalysis
{
	seq_t sequence;
	bifa_hits_t results;

	typedef boost::shared_ptr< BiFaAnalysis > ptr_t;
	typedef std::map< std::string, ptr_t > map_t;

protected:
	friend class boost::serialization::access;
    template< typename Archive >
    void serialize( Archive & ar, const unsigned int version )
    {
        ar & results;
		ar & sequence;
    }
};

void bifa_hits_2_match_results(
	const bifa_hits_t & hits,
	match_result_vec_t & results,
	double threshold = 0.0 );


BIO_NS_END

#endif //BIO_BIFA_ANALYSIS_H_

