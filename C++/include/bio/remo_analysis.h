#ifndef BIO_REMO_ANALYSIS_H_
#define BIO_REMO_ANALYSIS_H_

#include "bio/defs.h"
#include "bio/remo.h"
#include "bio/run_match.h"
#include "bio/binding_model.h"
#include "bio/bifa_analysis.h"


#include <sstream>


namespace boost { namespace program_options {
	//forward decl
	class options_description;
} }



BIO_NS_START




struct ReMoAnalysis
{
	typedef boost::shared_ptr< ReMoAnalysis > ptr_t;
	typedef std::map< ReMoRange, ptr_t > sequence_map_t; /**< A map that holds the analysis for a given sequence. */
	typedef std::map< ReMoLocation, sequence_map_t > location_map_t; /**< A map that holds the analysis for a given location in a sequence. */
	typedef std::map< std::string, location_map_t > map_t; /**< A map that holds the analysis for a set of sequences. */

	typedef boost::tuple< std::string, ReMoLocation, ReMoRange > range_key_t;

	bifa_hits_t results;
	seq_t sequence;

	/** Converts a remo map to a bifa map. */
	static void convert(map_t & remo_map, BiFaAnalysis::map_t & bifa_map);

protected:
	friend class boost::serialization::access;
    template<class Archive>
    void save( Archive & ar, const unsigned int version ) const
    {
		serialise< bifa_hits_t::value_type::binder_t >( ar, results );
		ar << sequence;
    }
    template<class Archive>
    void load( Archive & ar, const unsigned int version )
    {
		deserialise< bifa_hits_t::value_type::binder_t >( ar, results );
		ar >> sequence;
    }
    BOOST_SERIALIZATION_SPLIT_MEMBER()
};



template <class Fn>
void
for_each_double_hit(
	const ReMoAnalysis::map_t & analysis_map,
	unsigned distance,
	Fn fn,
	bool include_twins = false)
{
	unsigned num_sequences = 0;
	unsigned num_remos = 0;
	for (typename ReMoAnalysis::map_t::const_iterator s = analysis_map.begin();
		analysis_map.end() != s;
		++s, ++num_sequences)
	{
		//for each location
		for (typename ReMoAnalysis::location_map_t::const_iterator l = s->second.begin();
			s->second.end() != l;
			++l)
		{
			//for each remo
			for (typename ReMoAnalysis::sequence_map_t::const_iterator r = l->second.begin();
				l->second.end() != r;
				++r, ++num_remos)
			{
				//for each hit
				for (typename match_result_vec_t::const_iterator h1 = r->second.get()->results.begin();
					!( ( typename match_result_vec_t::const_iterator )( r->second.get()->results.end() ) == h1 );
					++h1)
				{
					const unsigned size_1 = BiobaseDb::singleton().get_pssm_entry(h1->link)->get_size();

					//look for other hits ahead of h1 but inside the distance
					for (typename match_result_vec_t::const_iterator h2 = r->second.get()->results.begin();
						!( ( typename match_result_vec_t::const_iterator )( r->second.get()->results.end() ) == h2 );
						++h2)
					{
						//check h2 is ahead of h1
						if (h2->result.position <= h1->result.position)
						{
							//it isn't
							continue;
						}

						//check h2 is not too far ahead, i.e. inside the distance
						if (unsigned(h2->result.position - h1->result.position) > distance + size_1)
						{
							//it is too far ahead
							continue;
						}

						//check they don't overlap
						if (unsigned(h2->result.position - h1->result.position) <= size_1)
						{
							continue;
						}

						//check they are not the same pssm
						if (! include_twins && h1->link == h2->link)
						{
							//they are
							continue;
						}

						fn(s->first, *h1, *h2);
					}
				}
			}
		}
	}

	std::cout << num_remos << " remos in " << num_sequences << " sequences\n";
}



class AnalysisVisitor
{
private:
	std::string analysis_filename;
	std::string sequence_name_regex;
	bool is_bifa_analysis;

	BiFaAnalysis::map_t bifa_analysis;
	ReMoAnalysis::map_t remo_analysis;

protected:
	unsigned num_sequences; /**< The number of sequences visited. */
	unsigned num_remos; /**< The number of remos visited. */
	unsigned num_bases; /**< The number of bases visited. */
	unsigned num_hits; /**< The number of putative binding sites. */
	double expected_num_hits; /**< The expected number of binding sites. */

public:
	void add_analysis_options(boost::program_options::options_description & options);

	/** Called before visiting the remos in each sequence group. Return false to ignore remos in this group. */
	virtual bool visit_sequence_group(const std::string & seq_group_name);

	/** Called after visiting the remos in each sequence group. */
	virtual void leave_sequence_group(const std::string & seq_group_name);

	/** Override to implement functionality for each remo visited. */
	virtual void visit_remo(
		const std::string & seq_group_name,
		ReMoLocation location,
		const ReMoRange & range,
		const std::string & remo_name,
		bifa_hits_t & results,
		const seq_t & sequence);

	void deserialise_analysis();
	void serialise_analysis(const std::string & filename) const;
	void visit_remo_analysis(bool show_progress = true);
};



BIO_NS_END

#endif //BIO_REMO_ANALYSIS_H_

