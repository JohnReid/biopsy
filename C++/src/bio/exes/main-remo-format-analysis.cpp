/* Copyright John Reid 2007
*/

#include "bio-pch.h"




#include "bio/application.h"
#include "bio/remo.h"
#include "bio/remo_analysis.h"
#include "bio/biobase_binding_model.h"
USING_BIO_NS

namespace po = boost::program_options;

#include <exception>
#include <iostream>
#include <fstream>
#include <string>
using namespace std;






struct FormatAnalysisApp : Application, AnalysisVisitor
{
	std::string output;
	std::ostream * output_stream;
	double threshold;

	FormatAnalysisApp()
	{
		get_options().add_options()
			( "output,o", po::value( &output )->default_value( "" ), "output filename" )
			( "threshold,t", po::value( &threshold )->default_value( 0.0 ), "threshold" )
			;

		add_analysis_options( get_options() );
	}

	bool visit_sequence_group(const std::string & seq_group_name)
	{
		*output_stream << "GENE," << seq_group_name << '\n';
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
		*output_stream
			<< "REMO,"
			<< location
			<< ","
			<< range.start
			<< ","
			<< range.end
			<< ","
			<< sequence
			<< "\n";

		typedef bifa_hits_t::index< BindingHitSet< BindingModel >::position >::type by_position_t;
		const by_position_t & by_position = results.get< BindingHitSet< BindingModel >::position >();
		for( by_position_t::const_iterator h = by_position.begin();
			by_position.end() != h;
			++h )
		{
			if( h->get_p_binding() >= threshold )
			{
				*output_stream 
					<< h->get_binder() 
					<< ","
					<< h->get_position() 
					<< ","
					<< ( h->is_complementary() ? "+" : "-" )
					<< ","
					<< h->get_length() 
					<< ","
					<< h->get_p_binding() 
					<< '\n';
			}
		}
	}

	int task()
	{
		boost::scoped_ptr< std::ofstream > output_file_stream;
		if( "" != output )
		{
			std::cout << "Writing output to \"" << output << "\"\n";
			output_file_stream.reset( new std::ofstream( output.c_str() ) );
			if( ! output_file_stream )
			{
				throw std::logic_error( BIO_MAKE_STRING( "Could not open \"" << output << "\"" ) );
			}
			output_stream = output_file_stream.get();
		}
		else
		{
			output_stream = &cout;
		}

		deserialise_analysis();
		visit_remo_analysis( false );

		return 0;
	}
};

int
main(int argc, char * argv[])
{
	return FormatAnalysisApp().main(argc, argv);
}

