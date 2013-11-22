/* Copyright John Reid 2007
*/

#include "bio-pch.h"




#include "bio/application.h"
#include "bio/remo.h"
#include "bio/remo_analysis.h"
#include "bio/svg_match.h"
#include "bio/file.h"
#include "bio/pathway_associations.h"
#include "bio/biobase_binding_model.h" //need this to register BIO_NS::BiobaseBindingModel::parameter_t class!
USING_BIO_NS

#include <boost/filesystem/path.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/operations.hpp>
using namespace boost;
namespace po = boost::program_options;
namespace fs = boost::filesystem;

#include <iostream>
#include <fstream>
#include <string>
using namespace std;




struct Analysis2SvgApp : Application, AnalysisVisitor
{
	BIO_NS::float_t threshold;
	bool open_svg;
	std::string prefix;

	Analysis2SvgApp()
		: open_svg( false )
	{
		get_options().add_options()
			("threshold,t", po::value(&threshold)->default_value(0.05f), "threshold")
			("open_svg,o", po::bool_switch(&open_svg), "open svg files")
			("prefix,p", po::value(&prefix)->default_value("seq"), "prefix")
			;

		add_analysis_options(get_options());
	}

	void visit_remo(
		const std::string & seq_group_name,
		ReMoLocation location,
		const ReMoRange & range,
		const std::string & remo_name,
		bifa_hits_t & hits,
		const seq_t & sequence)
	{
		match_result_vec_t results;
		bifa_hits_2_match_results( hits, results, threshold );
		if( results.empty() )
		{
			return;
		}

		//cout << "Creating SVG for analysis of " << remo_name << "\n";

		//make sure we have an SVG directory
		fs::path svg_dir( "SVG" );
		if( ! fs::exists( svg_dir ) )
		{
			fs::create_directory( svg_dir );
		}

		//get a filename that does not exist already
		std::string filename;
		unsigned index = 0;
		do
		{
			filename = BIO_MAKE_STRING( prefix << "_" << remo_name << "_" << index << ".svg" );
			BIO_FILENAMEIFY(filename);
			++index;
		}
		while (	fs::exists( svg_dir / fs::path( filename ) ) );

		//build and show the svg
		build_svg(
			svg_dir / filename,
			remo_name,
			sequence,
			threshold,
			results,
			16,
			false,
			open_svg );
	}

	int task()
	{
		deserialise_analysis();

		//load biobase and pathway info
		BiobaseDb::singleton();
		PathwayAssociations::singleton();

		visit_remo_analysis();

		return 0;
	}
};

int
main(int argc, char * argv[])
{
	return Analysis2SvgApp().main(argc, argv);
}

