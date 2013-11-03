/* Copyright John Reid 2007
*/

#include "bio-pch.h"




#include "bio/application.h"
#include "bio/remo_analysis.h"
USING_BIO_NS

namespace po = boost::program_options;

#include <iostream>
#include <fstream>
using namespace std;



struct AnalysisFilterApp : Application, AnalysisVisitor
{
	typedef std::vector< std::string > string_vec_t;

	std::string output_analysis_filename;

	AnalysisFilterApp()
	{
		get_options().add_options()
			("output,o", po::value(&output_analysis_filename), "output analysis file")
			;

		add_analysis_options(get_options());

		get_positional_options().add("motif", -1);
	}

	void visit_remo(
		const std::string & seq_group_name,
		ReMoLocation location,
		const ReMoRange & range,
		const std::string & remo_name,
		match_result_vec_t & results,
		const seq_t & sequence)
	{
	}




	int task()
	{
		deserialise_analysis();

		visit_remo_analysis();

		//do we want to write the output?
		if ("" != output_analysis_filename)
		{
			//serialise
			serialise_analysis(output_analysis_filename);
		}

		return 0;
	}
};

int
main(int argc, char * argv[])
{
	return AnalysisFilterApp().main(argc, argv);
}

