/* Copyright John Reid 2007
*/

#include "bio-pch.h"




#include "bio/application.h"
#include "bio/counter.h"
#include "bio/remo.h"
#include "bio/remo_analysis.h"
#include "bio/biobase_binding_model.h"
USING_BIO_NS

#include <boost/filesystem/path.hpp>
#include <boost/filesystem/fstream.hpp>
using namespace boost;
namespace po = boost::program_options;
namespace fs = boost::filesystem;

#include <iostream>
#include <fstream>
#include <string>
using namespace std;



struct AnalysisInfoApp : Application, AnalysisVisitor
{
	typedef Counter< unsigned > p_binding_counter_t;

	unsigned num_sequences;
	unsigned num_remos;
	unsigned num_hits;
	unsigned num_bases;
	p_binding_counter_t p_binding_counter;

	AnalysisInfoApp()
	{
		add_analysis_options(get_options());
	}

	bool visit_sequence_group(const std::string & seq_group_name)
	{
		++num_sequences;
		return true;
	}

	void visit_remo(
		const std::string & seq_group_name,
		ReMoLocation location,
		const ReMoRange & range,
		const std::string & remo_name,
		bifa_hits_t & results,
		const seq_t & sequence)
	{
		++num_remos;
		num_hits += results.size();
		num_bases += sequence.size();

		BOOST_FOREACH( const bifa_hits_t::value_type & hit, results )
		{
			const unsigned p_binding = unsigned( 100 * hit.get_p_binding() );
			p_binding_counter.increment( p_binding );
		}
	}

	int task()
	{
		deserialise_analysis();

		//for each sequence...
		num_sequences = 0;
		num_remos = 0;
		num_hits = 0;
		num_bases = 0;

		visit_remo_analysis( true );

		p_binding_counter.print(
			false,
			std::cout,
			"p(binding)",
			true,
			0,
			60,
			true );

		cout << "# genes: " << num_sequences << "\n";
		cout << "# remos: " << num_remos << "\n";
		cout << "# hits: " << num_hits << "\n";
		cout << "# bases: " << num_bases << "\n";

		return 0;
	}
};

int
main(int argc, char * argv[])
{
	return AnalysisInfoApp().main(argc, argv);
}


